#include "DualParticleSphLiquid.h"

#include "Geometries/ImplicitSurface.h"
#include "Utilities/Constants.h"

namespace PhysX {

    template<int Dim> void DualParticleSphLiquid<Dim>::writeDescription(YAML::Node & root) const {
        { // Description of particles.
            YAML::Node node;
            node["name"]                       = "particles";
            node["data_mode"]                  = "dynamic";
            node["primitive_type"]             = "point_list";
            node["material"]["diffuse_albedo"] = (Vector4f(52, 108, 156, 255) / 255).eval(); // Haijun Blue
            node["indexed"]                    = false;
            node["color_map"]["enabled"]       = true;
            root["objects"].push_back(node);
        }
        { // Description of particles.
            YAML::Node node;
            node["name"]                       = "virtual_particles";
            node["data_mode"]                  = "dynamic";
            node["primitive_type"]             = "point_list";
            node["material"]["diffuse_albedo"] = (Vector4f(52, 108, 156, 255) / 255).eval(); // Haijun Blue
            node["indexed"]                    = false;
            node["color_map"]["enabled"]       = true;
            root["objects"].push_back(node);
        }
        { // Description of particles.
            YAML::Node node;
            node["name"]                       = "boundary_particles";
            node["data_mode"]                  = "dynamic";
            node["primitive_type"]             = "point_list";
            node["material"]["diffuse_albedo"] = (Vector4f(52, 108, 156, 255) / 255).eval(); // Haijun Blue
            node["indexed"]                    = false;
            node["color_map"]["enabled"]       = false;
            root["objects"].push_back(node);
        }
    }

    template<int Dim>
    void DualParticleSphLiquid<Dim>::writeFrame(const std::string & frameDir, const bool staticDraw) const {
        { // Write particles.
            std::ofstream fout(frameDir + "/particles.mesh", std::ios::binary);
            IO::writeValue(fout, uint(_particles.size()));
            _particles.forEach(
                [&](const int i) { IO::writeValue(fout, _particles.positions[i].template cast<float>().eval()); });
            if constexpr (Dim == 3) {
                _particles.forEach([&](const int i) { IO::writeValue(fout, VectorDf::Unit(2).eval()); });
            }
            _particles.forEach([&](const int i) { IO::writeValue(fout, float(_velocities[i].norm())); });
        }
        { // Write particles.
            std::ofstream fout(frameDir + "/virtual_particles.mesh", std::ios::binary);
            IO::writeValue(fout, uint(_virtual_particles.size()));
            _virtual_particles.forEach([&](const int i) {
                IO::writeValue(fout, _virtual_particles.positions[i].template cast<float>().eval());
            });
            if constexpr (Dim == 3) {
                _virtual_particles.forEach([&](const int i) { IO::writeValue(fout, VectorDf::Unit(2).eval()); });
            }
            _virtual_particles.forEach(
                [&](const int i) { IO::writeValue(fout, float(_virtual_particles.volumes[i])); });
        }
        { // Write particles.
            std::ofstream fout(frameDir + "/boundary_particles.mesh", std::ios::binary);
            IO::writeValue(fout, uint(_boundary_particles.size()));
            _boundary_particles.forEach([&](const int i) {
                IO::writeValue(fout, _boundary_particles.positions[i].template cast<float>().eval());
            });
            if constexpr (Dim == 3) {
                _boundary_particles.forEach([&](const int i) { IO::writeValue(fout, VectorDf::Unit(2).eval()); });
            }
            //_particles.forEach([&](const int i) { IO::writeValue(fout, float(_velocities[i].norm())); });
        }
    }

    template<int Dim> void DualParticleSphLiquid<Dim>::initialize() {
        _boundary_particles.addSurface(
            ComplementarySurface<Dim>(std::make_unique<ImplicitBox<Dim>>(-VectorDr::Ones(), 2 * VectorDr::Ones())),
            VectorDr::Zero(),
            VectorDr::Ones() * 2 * 1.2);

        _boundary_velocity.resize(&_boundary_particles);

        _boundary_particles.setMass(_particles.mass());

        _particles.resetNearbySearcher();
        _particles.computeInfo();

        _targetDensity = _particles.densities.max();

        _boundary_particles.resetNearbySearcher();
        _boundary_particles.computeInfo();

        _virtual_particles.setAlpha(_particles, _targetDensity);

        _virtual_particles.setKappa(1. / _particles.getPackedKernelSum());

        /*double omega = 2.;

        if constexpr (Dim == 2)
            _particles.parallelForEach([&](const int i) {
                _velocities[i] = omega * _particles.positions[i].y() * VectorDr::Unit(0)
                    - omega * _particles.positions[i].x() * VectorDr::Unit(1);
            });*/

        std::cout << "Boundary Particle size: " << _boundary_particles.size() << std::endl;
        std::cout << "Virtual Particle size: " << _virtual_particles.size() << std::endl;
    }

    template<int Dim> void DualParticleSphLiquid<Dim>::advance(const real dt) {
        moveParticles(dt);

        applyExternalForces(dt);
        applyViscosityForce(dt);

        _virtual_particles.generateParticles(_particles);

        applyPressureForce(dt);
    }

    template<int Dim> void DualParticleSphLiquid<Dim>::applyPressureForce(const real dt) {
        _virtual_particles.project(
            _velocities, _boundary_velocity, _particles, _boundary_particles, _targetDensity, dt);
    }

    template class DualParticleSphLiquid<2>;
    template class DualParticleSphLiquid<3>;

} // namespace PhysX
