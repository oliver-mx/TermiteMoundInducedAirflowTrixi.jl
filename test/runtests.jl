using Trixi
using TermiteMoundInducedAirflowTrixi
using Test

# We run tests with CI jobs setting the `TRIXI_TEST` environment
# variable to determine the subset of tests to execute.
const TRIXI_TEST = get(ENV, "TRIXI_TEST", "all")

@time @testset "TermiteMoundInducedAirflowTrixi.jl tests" begin
    @time if TRIXI_TEST == "all" || TRIXI_TEST == "termite_mound"
        include("test_termite_mound.jl")
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "passive_house"
        include("test_passive_house.jl")
    end

    @time if TRIXI_TEST == "all" || TRIXI_TEST == "upstream"
        @testset "Namespace conflicts" begin
            # Test for namespace conflicts between TermiteMoundInducedAirflowTrixi.jl and Trixi.jl
            for name in names(Trixi)
                @test !(name in names(TermiteMoundInducedAirflowTrixi))
            end
        end
    
        # Run upstream tests for each mesh and dimension to test compatibility with Trixi.jl
        include("test_upstream.jl")
    end

end