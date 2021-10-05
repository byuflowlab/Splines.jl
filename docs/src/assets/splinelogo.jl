using Plots
import Splines

#circle
cps = [1  0 1/8 1;
       sqrt(2)/2  sqrt(2)/2 2*sqrt(2)/2/8 sqrt(2)/2;
       0  1 3/8 1;
      -sqrt(2)/2  sqrt(2)/2 4*sqrt(2)/2/8 sqrt(2)/2;
      -1  0 5/8 1;
      -sqrt(2)/2 -sqrt(2)/2 6*sqrt(2)/2/8 sqrt(2)/2;
       0 -1 7/8 1;
       sqrt(2)/2 -sqrt(2)/2 8*sqrt(2)/2/8 sqrt(2)/2;
       1  0 9/8 1]

knots = [0 0 0 1/4 1/4 1/2 1/2 3/4 3/4 1 1 1]
u = collect(0:0.01:1.0) #parametric points
p = 2 #curve order
n = length(cps[:, 1])-1
Cw = zeros(length(u), length(cps[1, :]))
for i = 1:length(u)
    Cw[i, :] = Splines.curvepoint(n, p, knots, cps, u[i])
end

plot3d(Cw[:, 1], Cw[:, 2], Cw[:,3], linewidths=10, aspectratio=:equal, grid=:off, label="")


