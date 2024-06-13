using Test

@testset "HarmonicVariable Tests" begin
    # Test display function
    @testset "Display" begin
        var = HarmonicVariable(:x, "x", "u", 1.0, 2.0)
        @test display(var) == :x
        vars = [
            HarmonicVariable(:x, "x", "u", 1.0, 2.0),
            HarmonicVariable(:y, "y", "v", 2.0, 3.0),
        ]
        @test display(vars) == [:x, :y]
    end

    # Test _coordinate_transform function
    @testset "_coordinate_transform" begin
        @test _coordinate_transform(2.0, 1.0, 0.0, "u") == 2.0
        @test _coordinate_transform(2.0, 1.0, 0.0, "v") == 0.0
        @test _coordinate_transform(2.0, 1.0, 0.0, "a") == 2.0
    end

    # Test _create_harmonic_variable function
    @testset "_create_harmonic_variable" begin
        rule, hvar = _create_harmonic_variable(2.0, 1.0, 0.0, "u"; new_symbol="x")
        @test rule == 2.0
        @test hvar == HarmonicVariable(:x, "u_{2.0}", "u", 1.0, 2.0)
    end

    # Test substitute_all function
    @testset "substitute_all" begin
        eq = 2.0 * :x + 3.0 * :y
        rules = Dict(:x => 1.0, :y => 2.0)
        @test substitute_all(eq, rules) == 2.0 + 6.0
        var = HarmonicVariable(:x, "x", "u", 1.0, 2.0)
        @test substitute_all(var, rules) == HarmonicVariable(1.0, "x", "u", 1.0, 2.0)
        vars = [
            HarmonicVariable(:x, "x", "u", 1.0, 2.0),
            HarmonicVariable(:y, "y", "v", 2.0, 3.0),
        ]
        @test substitute_all(vars, rules) == [
            HarmonicVariable(1.0, "x", "u", 1.0, 2.0),
            HarmonicVariable(2.0, "y", "v", 2.0, 3.0),
        ]
    end

    # Test Symbolics.get_variables function
    @testset "Symbolics.get_variables" begin
        vars = [1.0, 2.0, 3.0]
        @test Symbolics.get_variables(vars) == [1.0, 2.0, 3.0]
        var = HarmonicVariable(:x, "x", "u", 1.0, 2.0)
        @test Symbolics.get_variables(var) == [:x]
    end

    # Test isequal function
    @testset "isequal" begin
        var1 = HarmonicVariable(:x, "x", "u", 1.0, 2.0)
        var2 = HarmonicVariable(:x, "x", "u", 1.0, 2.0)
        var3 = HarmonicVariable(:y, "y", "v", 2.0, 3.0)
        @test isequal(var1, var2) == true
        @test isequal(var1, var3) == false
    end
end
