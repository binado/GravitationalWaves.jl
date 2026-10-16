@testset "model building blocks" begin
    G = GWWaveforms

    @testset "TaylorF2 phasing" begin
        eta = 36.0 * 29.0 / 65.0^2
        pn = G.taylorf2_pn_phasing(36.0, 29.0, 0.0, 0.0)
        @test pn.v[1] ≈ 3 / (128 * eta)
        @test pn.v[2] == 0
        # Without spin terms the spins do not matter.
        @test G.taylorf2_pn_phasing(36.0, 29.0, 0.5, -0.3; spin_order = 0) == pn
        @test G.taylorf2_pn_phasing(36.0, 29.0, 0.5, -0.3) != pn
        # Truncation to the Newtonian order keeps only the leading term.
        newtonian = G.truncate_phasing(G.taylorf2_pn_phasing(36.0, 29.0, 0.5, 0.3), 0)
        @test G.taylorf2_phasing(newtonian, 0.2) ≈ newtonian.v[1] / 0.2^5
        # The higher orders make a small, positive correction to the Newtonian phasing.
        full = G.taylorf2_pn_phasing(36.0, 29.0, 0.0, 0.0)
        @test 0 < abs(G.taylorf2_phasing(full, 0.18) / (full.v[1] / 0.18^5) - 1) < 0.2
    end
end
