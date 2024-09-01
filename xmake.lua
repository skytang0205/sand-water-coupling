set_project("VCL-Physx")
set_version("0.0")
set_xmakever("2.3.9")

-- Usage:
--   <plain>
--   > xmake && xmake b examples && xmake r LevelSetLiquidTest
--   > xmake r vcl-viewer
--   <openmp>
--   > xmake f --openmp=true && xmake && xmake b examples && xmake r LevelSetLiquidTest
--   > xmake r vcl-viewer
--   optional: install to some directory (such as ./dist)
--   > xmake install -o dist && xmake install -o dist examples

add_requires("eigen", "glfw", "glad", "yaml-cpp", "fmt")
pkgs = {"eigen", "glfw", "glad", "yaml-cpp", "fmt"}

add_rules("mode.debug", "mode.release")
add_cxflags("/utf-8")

-- BLAS is used for only dense matrix operations, which is not useful for this project.
-- option("use_blas")
--     set_default(false)
--     set_showmenu(true)
--     set_description("Enable blas support.")
--     add_defines("EIGEN_USE_BLAS")
-- option_end()

-- if has_config("use_blas") then
--     add_requires("openblas")
-- end

option("openmp")
    set_default(false)
    set_showmenu(true)
    set_description("Enable openmp support.")
    if is_plat("windows") then
        add_cxflags("/openmp")
    else
        add_cxflags("-fopenmp")
    end
    add_defines("_OPENMP")
    add_defines("EIGEN_HAS_OPENMP")
-- benchmark: plain ~ 132.72s, openmp ~ 178.86s
option_end()

option("use_float")
    set_default(false)
    set_showmenu(true)
    set_description("Use single-precision floating point numbers instead of double-precision ones.")
    add_defines("USE_FLOAT")
option_end()

option("common")
    set_default(true)
    set_showmenu(false)
    set_languages("cxx20")
    if is_plat("windows") then
        add_cxxflags("/MT")
    end
option_end()

target("vcl-physx")
    set_default(true)
    set_kind("$(kind)")
    add_options("common")
    add_options("openmp")
    -- add_options("use_blas")
    add_packages(unpack(pkgs))
    add_includedirs(path.join(os.curdir(), "Cores"), {public = true})
    local modules = {"Geometries", "Graphics", "Physics", "Solvers", "Structures", "Utilities"}
    for _, module in ipairs(modules) do
        add_files("Cores/" .. module .. "/*.cpp")
        add_headerfiles("Cores/(" .. module .. "/*.h)")
    end
target_end()

target("vcl-viewer")
    set_default(true)
    set_kind("binary")
    add_options("common")
    add_deps("vcl-physx")
    add_packages(unpack(pkgs))
    add_files("Cores/Viewer/*.cpp")
target_end()

local examples = {"EulerianFluidTest", "LevelSetLiquidTest", "ParticleInCellLiquidTest", "MatPointSubstancesTest", "SpringMassSystemTest", "SmthPartHydrodLiquidTest", "DEMParticleSandTest"}
for _, example in ipairs(examples) do

target(example)
    set_default(false)
    set_kind("binary")
    add_options("common")
    -- add_options("use_blas")
    add_packages(unpack(pkgs))
    add_deps("vcl-physx")
    add_files(path.join("Demos", example, "*.cpp"))
target_end()

end

target("examples")
    set_default(false)
    set_kind("phony")
    add_deps(unpack(examples))
target_end()

includes("Develop")
