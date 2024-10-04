set_project("fluid-grid-2d")
set_languages("cxxlatest")

add_rules("mode.debug", "mode.release")

includes("SimViewer")

add_requires("amgcl")
add_requires("cxxopts v3.1.1")
add_requires("eigen")
add_requires("spdlog", { configs = { fmt_external = true } })
add_requires("tbb")
add_requires("yaml-cpp 0.7.0")

target("demo")
    set_kind("binary")
    add_packages("amgcl")
    add_packages("cxxopts")
    add_packages("eigen" )
    add_packages("spdlog")
    add_packages("tbb")
    add_packages("yaml-cpp")

    add_includedirs("core")
    add_headerfiles("core/**.h")
    add_files("core/**.cpp")
    add_headerfiles("demo/**.h")
    add_files("demo/**.cpp")
