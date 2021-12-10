import HarmonicBalance: drop_powers, max_power
import HarmonicBalance.Symbolics.expand

@variables a,b,c 

@test max_power(a^2 + b, a) == 2
@test max_power( a*((a+b)^4 )^2 + a, a) == 9


@test isequal(drop_powers(a^2+b, a, 1) , b)
@test isequal(drop_powers( (a+b)^2, a,1), b^2)
@test isequal(drop_powers( (a+b)^2, [a,b],1), 0)

@test isequal(drop_powers((a+b)^3 + (a+b)^5 , [a,b],4), expand((a+b)^3))
