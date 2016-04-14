clear all
load('data_rgdp_and_volatility')

for n = 1:25
    figure(n)
    plot(deflator(:, n))
end