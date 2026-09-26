using Radiant, Test

@testset "2D first-order closure positivity and balance" begin
    mu = inv(sqrt(3.0))
    dd = reshape([-1.0,2.0],2,1,1)
    dg = reshape([0.0,1.0],2,1,1)
    for (dx,dy,weights,expected_cell,expected_y) in (
        (1e-4,1e-3,dd,1/11,-9/11),
        (1e-4,1e-4,dd,1/2,0.0),
        (1e-4,1e-3,dg,1/11,1/11),
    )
        cell,xout,yout,_ = Radiant.flux_2D_BFP(mu,mu,0.0,0.0,0.0,[0.0],4.0,dx,dy,
            [0.0],[0.0],[1.0],[0.0],1,1,1,[1.0],dg,weights,weights,false,zeros(1,1,1),true)
        @test cell[1] ≈ expected_cell
        @test yout[1] ≈ expected_y atol=1e-14
        # Zero volume source and removal: outgoing minus incoming balances.
        @test (mu/dx)*xout[1] + (mu/dy)*(yout[1]-1.0) ≈ 0.0 atol=1e-10
    end
end
