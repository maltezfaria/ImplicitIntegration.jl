using SafeTestsets

@safetestset "Aqua tests (code quality)" begin
    include("aqua_test.jl")
end

@safetestset "Subfunction tests" begin
    include("subfunction_test.jl")
end

@safetestset "Integration tests" begin
    include("integration_test.jl")
end