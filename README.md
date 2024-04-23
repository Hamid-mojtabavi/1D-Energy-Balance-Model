# 1D Energy Balance Model
## Introduction

1 - dimensional energy balance model based on the text book 
by Henderson - Sellers and McGuffie (1987):

Implemented in Python by Rafaila Grigoriou, rafaila.gr@gmail.com***
Updated by Dr. Seyedhamidreza Mojtabavi (Seyedhamidreza.Mojtabavi@gmail.com)

Two - hemispheres version (18 latitude zones)

This program computes the surface temperature of a latitude zone from the local
energy balance.It includes the ice-albedo feedback

- reads the length of the run and the fraction of the solar constant from console
- reads any additional radiative forcings from file
- writes the final temperature and energy flux distributions to console as well as a
snapshot file
- writes a time series of global mean temperature to a history file

The code includes a function for the identification of the value of fractional solar constant just for completely glaciating the Earth (ice edge at 0^o N) starting from an unglaciated Earth and vice versa (value of fractional solar constant just for unglaciating Earth, starting from a completely glaciated Earth)

## Running examples

Calling the EBM function for the default values of the input parameters:

- scon0 = solar constant (W m^(-2)) (optional, def: 1370.0)
- tzero = pure water freezing point in K (optional, def: 273.15)
- tcrit = temperature at which the surface becomes ice covered (optional, def: -10)
- t0 = initial temperature (deg C)(optional, def: 15.0)
- sfrac = fractional solar constant (optional, def: 1.0)
- ktrans = transport coefficient (W m^(-2) K^(-1)) -  relaxation coeffient in Newtonian (linear) approximation of meridional heat transport divergence (optional, def: 3.81)
- aice = albedo of snow and ice (optional, def: 0.62)
- alw = constant part of linearized longwave radiation (W m^(-2))(optional, def: 204.0)
- blw = linear part of linearized longwave radiation (W m^(-2) K^(-1)) (optional, def: 2.17)
- to_print = boolean, when True the function prints model output in console (optional, def: False)


```python 
  df_EBM, tmean, yice_s, yice_n, jice_s, jice_n = 
  EnergyBalanceModel1D(ittmax, grid_step, asfc, year0, year_len, dt, cocean, s, latzone, to_print=True) 
```

Finding "Glaciated Earth Edge Point" and plot graph of temperature vs sfrac, going from non-glaciated Earth to glaciated Earth

```python
    sfrac_range = [0.5, 1.8]
    step = 0.01

    sfrac_ice, df_sfrac = find_sfrac_ice(sfrac_range, step, ittmax, grid_step, asfc, year0, year_len, dt, cocean, s, latzone,
                                         scon0, tzero, tcrit,  t0, ktrans, aice, alw, blw)
    print("Fractional solar constant value at Glaciated Earth Edge Point is:", sfrac_ice)
```

The main includes more running examples.

