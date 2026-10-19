using Test, GWSEOBNR, Enzyme
const S = GWSEOBNR

function validated_loss(x, retcode)
    S.check_evolution_retcode(retcode)
    return sum(abs2, x)
end

function checked_gradient(retcode)
    x = [0.4, -0.7]
    dx = zeros(2)
    Enzyme.autodiff(
        Enzyme.Reverse, Enzyme.Const(validated_loss), Enzyme.Active,
        Enzyme.Duplicated(x, dx), Enzyme.Const(retcode)
    )
    return dx
end

@testset "inactive solver-status validation" begin
    @test checked_gradient(S.SciMLBase.ReturnCode.Terminated) ≈ [0.8, -1.4]
    # The inactive helper must execute and reject invalid solves under AD too.
    @test_throws "EOB evolution exhausted time coverage without a physical guard" checked_gradient(S.SciMLBase.ReturnCode.Success)
    @test_throws "EOB evolution failed:" checked_gradient(S.SciMLBase.ReturnCode.MaxIters)
end
