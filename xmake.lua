set_project("VCL-Physx")
set_version("0.0")
set_xmakever("2.3.9")

add_requires("eigen", "glfw", "glad", "yaml-cpp", "fmt")
local pkgs = {"eigen", "glfw", "glad", "yaml-cpp", "fmt"}

add_rules("mode.debug", "mode.release")

--[[option("use_blas")
    set_default(false)
    set_showmenu(true)
    set_description("Enable blas support.")
    add_defines("EIGEN_USE_BLAS")
option_end()

if has_config("use_blas") then
    add_requires("openblas")
else
    add_requires("openblas", {optional = true})
end]]--

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

option("use_float")
    set_default(false)
    set_showmenu(true)
    set_description("Use single-precision floating point numbers instead of double-precision ones.")
    add_defines("USE_FLOAT")

target("vcl-physx")
    set_default(true)
    set_kind("$(kind)")
    set_languages("cxx20")
    if is_plat("windows") then
        add_cxxflags("/MT")
    end
    add_options("openmp")
    add_options("use_blas")
    add_packages(unpack(pkgs))
    --add_packages("openblas", {optional = true})
    add_includedirs(path.join(os.curdir(), "Cores"), {public = true})
    local modules = {"Geometries", "Graphics", "Physics", "Solvers", "Structures", "Utilities"}
    for _, module in ipairs(modules) do
        add_files("Cores/" .. module .. "/*.cpp")
        add_headerfiles("Cores/(" .. module .. "/*.h)")
    end

target("vcl-viewer")
    set_default(true)
    set_kind("binary")
    set_languages("cxx20")
    if is_plat("windows") then
        add_cxxflags("/MT")
    end
    add_deps("vcl-physx")
    add_packages(unpack(pkgs))
    add_files("Cores/Viewer/*.cpp")
    on_run(function (target)
        os.execv(target:targetfile(), {"-o", os.curdir() .. "/output"})
    end)

local examples = {"EulerianFluidTest", "LevelSetLiquidTest", "ParticleInCellLiquidTest", "SpringMassSystemTest"}
for _, example in ipairs(examples) do
    target(example)
        set_default(false)
        set_kind("binary")
        set_languages("cxx20")
        if is_plat("windows") then
            add_cxxflags("/MT")
        end
        --add_options("use_blas")
        add_packages(unpack(pkgs))
        --add_packages("openblas", {optional = true})
        add_deps("vcl-physx")
        add_files(path.join("Demos", example, "*.cpp"))
        on_run(function (target)
            os.execv(target:targetfile(), {"-o", os.curdir() .. "/output"})
        end)
end

target("examples")
    set_default(false)
    set_kind("phony")
    add_deps(unpack(examples))

task("view")
    on_run(function (name)
        print(os.curdir())
    end)
