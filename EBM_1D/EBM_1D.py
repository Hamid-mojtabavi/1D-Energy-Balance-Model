import pandas as pd
import mpmath as mp
from mpmath import mpf
import matplotlib.pyplot as plt
import numpy as np




mp.dps = 64  # increasing accuracy

'''
1 - dimensional energy balance model based on the text book 
by Henderson - Sellers and McGuffie (1987):

Two - hemispheres version (18 latitude zones)

This program computes the surface temperature of a latitude zone from the local
energy balance.It includes the ice-albedo feedback

- reads the length of the run and the fraction of the solar constant from console
- reads any additional radiative forcings from file
- writes the final temperature and energy flux distributions to console as well as a
snapshot file
- writes a time series of global mean temperature to a history file

***Implemented in Python by Rafaila Grigoriou, rafaila.gr@gmail.com***
'''

# DEFINITION OF CONSTANTS
pi = mpf('3.14159265358979')
deg2rad = pi / mpf(180)  # Degrees to radians conversion factor
deg2dist = mpf('111194.9')  # Distance that corresponds to one degree of latitude in units of m

jmt = 20  # number of grid points in the latitudinal direction
jmtm1 = jmt - 1


def albedo(Yt, Yu, tcrit, asfc, aice, t, apln):

    ''' This function computes the planetary albedo of the latitude zones.

    Input arguments:
        yt  = T grid cell latitude (deg N)
        yu  = U grid cell latitude (deg N)
        tcrit  = temperature at which the surface becomes ice covered
        asfc = surface albedo
        aice = albedo of snow and ice
        t  = temperature (deg C)
        apln = planetary albedo

    Output:
        apln = planetary albedo (actual value)
        yice_s = latitude of ice boundary (southern hemisphere)
        yice_n = latitude of ice boundary (northern hemisphere)
        jice_s = latitudinal index of ice boundary (southern hemisphere)
        jice_n = latitudinal index of ice boundary (northern hemisphere)

    References: Henderson-Sellers and McGuffie (1987)'''

    # dyt = Yt['dyt']
    # dyu = Yu['dyu']
    yt = Yt['yt']
    yu = Yu['yu']

    half_jmt = int(jmt/2)

    # ----------------- SOUTHERN HEMISPHERE

    # Find latitude of ice boundary

    for j in range(1, half_jmt):  # starting at South Pole
        if (t[j] < tcrit) & (t[j+1] >= tcrit):  # Actual grid cell colder than critical temperature, but grid cell to the north not: partial ice cover (assumes temperature monotonous function of latitude)
            # print("SOUTHERN: Actual grid cell colder than critical temperature, but grid cell to the north not: partial ice cover (assumes temperature monotonous function of latitude)")
            dphi = (-tcrit + t[j + 1]) / (t[j + 1] - t[j]) * 0.1745 / deg2rad
            yice_s = yt[j + 1] - dphi

            if yice_s > yu[j]:  # Ice boundary less than half - way between grid points
                # print("SOUTHERN: Ice boundary less than half - way between grid points")
                ysouth = yu[j]
                ynorth = yu[j + 1]
                wice = (mp.sin(yice_s * deg2rad) - mp.sin(ysouth * deg2rad)) / (mp.sin(ynorth * deg2rad) - mp.sin(ysouth * deg2rad))
                amix_s = asfc[j + 1] * (mpf(1) - wice) + aice * wice
                jice_s = j + 2
            else:  # Ice boundary more than half - way between grid points
                # print("SOUTHERN: Ice boundary more than half - way between grid points")
                ysouth = yu[j - 1]
                ynorth = yu[j]
                wice = (mp.sin(yice_s * deg2rad) - mp.sin(ysouth * deg2rad)) / (mp.sin(ynorth * deg2rad) - mp.sin(ysouth * deg2rad))
                amix_s = asfc[j] * (mpf(1) - wice) + aice * wice
                jice_s = j + 1

    # Set planetary albedo of zones
    if t[1] > tcrit:  # Totally ice - free case
        # print("SOUTHERN: Totally ice - free case")
        jice_s = 0
        yice_s = mpf(-90)

        for j in range(1, half_jmt):
            apln[j] = asfc[j]

    elif t[half_jmt-1] <= tcrit:  # Totally ice - covered case
        # print("SOUTHERN: Totally ice - covered case")
        jice_s = half_jmt
        yice_s = mpf(0)

        for j in range(1, half_jmt):
            apln[j] = aice

    else:  # Partially ice - covered case
        # print("SOUTHERN: Partially ice - covered case")
        for j in range(1, jice_s-1):
            apln[j] = aice

        apln[jice_s-1] = amix_s

        for j in range(jice_s, half_jmt):
            apln[j] = asfc[j]

    # ----------------- NORTHERN HEMISPHERE

    # Find latitude of ice boundary
    for j in range(jmtm1-1, half_jmt-1, -1):  # starting at North Pole
        if (t[j] < tcrit) & (t[j-1] >= tcrit): # Actual grid cell colder than critical temperature, but grid cell to the south not: partial ice cover (assumes temperature monotonous function of latitude)
            # print("NOTHERN: Actual grid cell colder than critical temperature, but grid cell to the south not: partial ice cover (assumes temperature monotonous function of latitude")
            dphi = (-tcrit + t[j-1])/(t[j-1] - t[j])*mpf('0.1745')/deg2rad
            yice_n = yt[j-1] + dphi

            if yice_n < yu[j-1]:  # Ice boundary less than half-way between grid points
               # print("NOTHERN: Ice boundary less than half-way between grid points")
               ysouth = yu[j-2]
               ynorth = yu[j-1]
               wice = (mp.sin(ynorth*deg2rad) - mp.sin(yice_n*deg2rad)) / (mp.sin(ynorth*deg2rad) - mp.sin(ysouth*deg2rad))
               amix_n = asfc[j-1]*(mpf(1) - wice) + aice*wice
               jice_n = j
            else:  # Ice boundary more than half-way between grid points
               # print("NOTHERN: Ice boundary more than half-way between grid points")
               ysouth = yu[j-1]
               ynorth = yu[j]
               wice = (mp.sin(ynorth*deg2rad) - mp.sin(yice_n*deg2rad)) / (mp.sin(ynorth*deg2rad) - mp.sin(ysouth*deg2rad))
               amix_n = asfc[j]*(mpf(1) - wice) + aice*wice
               jice_n = j+1

    # Set planetary albedo of zones
    if t[jmtm1-1] > tcrit:  # Totally ice-free case
        # print("NOTHERN: Totally ice-free case")
        jice_n = 0.0
        yice_n = mpf(90)

        for j in range(jmtm1-1, half_jmt-1, -1):
            apln[j] = asfc[j]

    elif t[half_jmt] <= tcrit:  # Totally ice-covered case
        # print("NOTHERN: Totally ice-covered case")
        jice_n = half_jmt
        yice_n = mpf(0)

        for j in range(jmtm1-1, half_jmt-1, -1):
            apln[j] = aice

    else:  # Partially ice-covered case
        # print("NOTHERN: Partially ice-covered case")
        for j in range(jmtm1-1, jice_n-1, -1):
            apln[j] = aice

        apln[jice_n-1] = amix_n

        for j in range(jice_n-2, half_jmt-1, -1):
            apln[j] = asfc[j]

    return apln, yice_s, yice_n, jice_s, jice_n


def EnergyBalanceModel1D(ittmax, grid_step, asfc, year0, year_len, dt, cocean, s, latzone,
                         scon0=1370.0, tzero=273.15, tcrit=-10.0, t0=15.0, sfrac=1.0, ktrans=3.81, aice=0.62,
                         alw=204.0, blw=2.17, to_print=False):

    ''' This function runs the 1D energy balance model, using the specified value of the input parameters.

        Input arguments:
            ittmax = maximum number of integer time counter (maximum number of timesteps)
            grid_step = grid step (10\B0)
            asfc = surface albedo
            year0 = initial model year
            year_len = number of days in year
            dt = time step in seconds
            tzero = pure water freezing point(K)
            cocean = heat capacity of the ocean (J m^(-2) K^(-1))
            s = solar distribution function
            latzone = latitude zones
            scon0 = solar constant (W m^(-2)) (optional, def: 1370.0)
            tzero = pure water freezing point in K (optional, def: 273.15)
            tcrit = temperature at which the surface becomes ice covered (optional, def: -10)
            t0 = initial temperature (deg C)(optional, def: 15.0)
            sfrac = fractional solar constant (optional, def: 1.0)
            ktrans = transport coefficient (W m^(-2) K^(-1)) -  relaxation coeffient in Newtonian (linear) approximation of meridional heat transport divergence (optional, def: 3.81)
            aice = albedo of snow and ice (optional, def: 0.62)
            alw = constant part of linearized longwave radiation (W m^(-2))(optional, def: 204.0)
            blw = linear part of linearized longwave radiation (W m^(-2) K^(-1)) (optional, def: 2.17)
            to_print = boolean, when True the function prints model output in console (optional, def: False)

        Output:
            df_ouput = dataframe with structure:
                'Zone': latzone (latitude zone)
                'Temp (deg C)': ti (Earth temperature per latitude zone)
                'Albedo': apln (planetary albedo - actual value)
                'LW (W m^(-2))': netlwpln (net longwave radiation absorbed by planet in W m^(-2))
                'SW (W m^(-2))': netswpln (net shortwave radiation absorbed by planet in W m^(-2))
            tmean = mean Earth temperature (deg C)
            yice_s = latitude of ice boundary (southern hemisphere)
            yice_n = latitude of ice boundary (northern hemisphere)
            jice_s = latitudinal index of ice boundary (southern hemisphere)
            jice_n = latitudinal index of ice boundary (northern hemisphere)

        *option2: output in dataframe (output_in_df = True)*
            Dataframe structure:



        References: Henderson-Sellers and McGuffie (1987) '''

    p2 = pi / mpf(2)  # "pi over 2"
    solin = sfrac * scon0 * 0.25 * s  # solar insolation (W m^(-2))

    # Set model grid:
    # Define T grid cells
    Yt = pd.DataFrame(columns=['yt', 'sint', 'cst'])
    Yu = pd.DataFrame(columns=['yu', 'sinu', 'csu'])

    for j in range(0, jmt):

        T_dict, U_dict = {}, {}
        T_dict['yt'] = mpf(-95)  # T grid cell latitude
        if j > 0:
            T_dict['yt'] = Yt['yt'][j-1] + mpf(10)
        T_dict['sint'] = mp.sin(T_dict['yt']*deg2rad)  # sine of T grid cell latitude
        T_dict['cst'] = mp.cos(T_dict['yt']*deg2rad)  # cosine of T grid cell latitude

        #Yt = Yt.append(T_dict, ignore_index=True)
        Yt = pd.concat([Yt, pd.DataFrame([T_dict])], ignore_index=True)

        # U_dict['dyu'] = mpf(10) * deg2dist  # U grid cell height
        U_dict['yu'] = mpf(-90)  # U grid cell latitude
        if j > 0:
            U_dict['yu'] = Yu['yu'][j-1] + mpf(10)

        U_dict['sinu'] = mp.sin(U_dict['yu'] * deg2rad)  # sine of U grid cell latitude
        U_dict['csu'] = mp.cos(U_dict['yu'] * deg2rad)  # cosine of U grid cell latitude

        #Yu = Yu.append(U_dict, ignore_index=True)
        Yu = pd.concat([Yu, pd.DataFrame([U_dict])], ignore_index=True)

    ti = [t0] * jmtm1  # Initialize temperature at old time step (deg C)
    tf = [None]*jmtm1  # Initialize temperature at old time step (deg C)
    apln = [None] * jmtm1

    for itt in range(0, ittmax):
        apln, yice_s, yice_n, jice_s, jice_n = albedo(Yt, Yu, tcrit, asfc, aice, ti, apln)  #  Compute planetary albedo of zones

        amean, tmean = 0.0, 0.0
        netswpln, netlwpln = [None]*jmtm1, [None]*jmtm1  # initialize net shortwave and longwave radiation absorbed by planet (W m^(-2))

        for j in range(1, jmtm1):
            a1 = p2 - (j - 1.) * grid_step  # northern boundary of grid cell
            a2 = a1 - grid_step  # southern boundary of grid cell
            amean = amean + (mp.sin(a1) - mp.sin(a2)) * 0.5 * apln[j]
            tmean = tmean + (mp.sin(a1) - mp.sin(a2)) * 0.5 * ti[j]

        # Update model time
        model_time = itt*dt/(86400.0*year_len) + year0

        for j in range(1, jmtm1):
            # Compute shortwave radiation budget
            dqtot = 0.0  # alternative to turn off radiative forcing from crowley file
            netswpln[j] = solin[j] * (1.0 - apln[j]) + dqtot

            # Compute longwave radiation budget
            tkelvin = ti[j] + tzero
            netlwpln[j] = alw + blw * ti[j]

            # Set source term for energy balance equation (W m^(-2))
            source = netswpln[j] - netlwpln[j]

            # Compute temperature at new time step
            tf[j] = ti[j] + (dt / cocean) * (-ktrans * (ti[j] - tmean) + source)

        ti = tf

    df_output = pd.DataFrame({'Zone': latzone, 'Temp (deg C)': ti, 'Albedo': apln, 'LW (W m^(-2))': netlwpln, 'SW (W m^(-2))': netswpln})
    df_output = df_output.iloc[1:, :]

    if to_print:
        print(df_output)
        print("The final global mean temperature is tmean =", tmean, " deg C")
        print("The latitude of the ice edge is at yice_s =", yice_s, "deg (jice_s =", jice_s, ")")
        print("                                   yice_n =", yice_n, "deg (jice_n =", jice_n, ")")

    return df_output, tmean, yice_s, yice_n, jice_s, jice_n


def find_sfrac_ice(sfrac_range, step, ittmax, grid_step, asfc, year0, year_len, dt, cocean, s, latzone,
                   scon0=1370.0, tzero=273.15, tcrit=-10.0, t0=15.0, ktrans=3.81, aice=0.62,
                   alw=204.0, blw=2.17, run_all=True):

    ''' This function runs an iterative loop of the EnergyBalanceModel1D for different values of the
    fractional solar constant on a defined range and step. The "Glaciated Earth Edge Point" (i.e. the value of sfrac
     needed to glaciate Earth completely) is identified. The mean temperature of the previous iteration is used as the
     initial temperature of the next iteration.

        Input arguments:
            sfrac_range = the range in which sfrac takes values (the direction of the range - i.e. from smaller
                          to larger or from larger to smaller sfrac values - defines the direction that the iterations
                          take place)
            step = the step of sfrac for the iterations
            EnergyBalanceModel1D input arguments
            run_all = boolean, when True all sfrac values in the range are used,
            when False the iterative algorithm stops once "Glaciated Earth Edge Point" is identified

        Output:
            sfrac_ice = the value of sfrac at the "Glaciated Earth Edge Point"
            df_sfrac = a dataframe with climate values for each sfrac value, with structure
                'sfrac': sfrac
                'tmean': tmean
                yice_n': yice_n
                'yice_s': yice_s           '''

    df_sfrac = pd.DataFrame(columns=['sfrac', 'tmean', 'yice_n', 'yice_s'])
    sfrac_ice = 0

    if sfrac_range[0] < sfrac_range[1]:  # sfrac values go from smaller to larger
        for sfrac in np.arange(sfrac_range[1], sfrac_range[0]-step, -step):
            df_EBM, tmean, yice_s, yice_n, jice_s, jice_n = \
                EnergyBalanceModel1D(ittmax, grid_step, asfc, year0, year_len, dt, cocean, s, latzone, scon0, tzero,
                                     tcrit, t0, sfrac, ktrans, aice, alw, blw, to_print=False)
            df_sfrac = df_sfrac.append({'sfrac': sfrac, 'tmean': tmean, 'yice_n': yice_n, 'yice_s': yice_s}, ignore_index=True)
            t0 = tmean
            if yice_n < 0.01:
                if sfrac_ice == 0:
                    sfrac_ice = sfrac
                if not run_all:
                    break
    else:  #sfrac values go from larger to smaller
        for sfrac in np.arange(sfrac_range[1], sfrac_range[0]+step, step):
            df_EBM, tmean, yice_s, yice_n, jice_s, jice_n = \
                EnergyBalanceModel1D(ittmax, grid_step, asfc, year0, year_len, dt, cocean, s, latzone, scon0, tzero,
                                     tcrit, t0, sfrac, ktrans, aice, alw, blw, to_print=False)
            df_sfrac = df_sfrac.append({'sfrac': sfrac, 'tmean': tmean, 'yice_n': yice_n, 'yice_s': yice_s}, ignore_index=True)
            t0 = tmean
            if yice_n > 0.0:
                if sfrac_ice == 0:
                    sfrac_ice = sfrac
                if not run_all:
                    break

    print("sfrac_ice:", sfrac_ice)
    df_sfrac = df_sfrac.astype(float)

    return sfrac_ice, df_sfrac


def plot_estimated_line(df, x, y, ax=None, vline_place=None, title='', xlabel='', ylabel='', grid=True, linewidth=1, color='blue', legend=True):
    ''' Support function used for plotting lines from dataframes, using markers for the given values and connect them with lines.
    The addition of a vertical dotted line is also possible.

        Input arguments:
            df = input pandas dataframe
            x = dependant variable name
            y = independant variable name
            ax = matplotlib.pyplot.axes object (optional, def: None)
            vline_place = when != None, defined the x location for the dotted vertical line (optional, def: None)
            title = graph title, string(optional, def: '')
            xlabel = graph xlabel, string (optional, def: '')
            ylabel = graph ylabel, string (optional, def: '')
            grid = boolean, when True adds a grid to the graph(optional, def: True)
            linewidth = the width of the line in the graph(optional, def: 1)
            color = the color of the line(optional, def: 'blue')
            legend = boolean, if True a legend is added to the graph (optional, def: True)

        Output:
            ax = matplotlib.pyplot.axes object            '''

    if ax is not None:
        df.plot(x=x, y=y, kind='scatter', title=title, xlabel=xlabel, ylabel=ylabel, grid=grid, ax=ax, color=color)
    else:
        ax = df.plot(x=x, y=y, kind='scatter', title=title, xlabel=xlabel, ylabel=ylabel, grid=grid, color=color)

    df.plot(x=x, y=y, kind='line', title=title, xlabel=xlabel, grid=grid, linewidth=linewidth, ax=ax, legend=legend, color=color)

    if vline_place is not None:
        plt.axvline(vline_place, color='grey', linewidth=2, linestyle='dashed')

    return ax


if __name__ == '__main__':

    # Set file names
    fnameinit = "init_2h.dat"

    df_init = pd.read_csv(fnameinit, skiprows=3, delimiter=r"\s+")
    for col in df_init.columns:
        df_init.rename({col: col.strip()})

    df_init.loc[-1] = [None, None, None, None]
    df_init.index = df_init.index + 1  # shifting index
    df_init.sort_index(inplace=True)

    # Set model parameters
    tzero = 273.15  # pure water freezing point (K)
    cocean = 2.94E08  # heat capacity of the ocean (J m^(-2) K^(-1)), with respect to top 70 m of the ocean (Hartmann 1994)
    scon0 = 1370.0  # solar constant (W m^(-2)), present - day value(Henderson - Sellers and McGuffie 1987)
    eccen = 0.0167  # numerical eccentricity of Earth's orbit, present - day value(Berger 1978)
    tcrit = -10.0  # temperature at which the surface becomes ice covered (def = -10.0)
    grid_step = mpf('3.14159') / mpf(18)  # grid step(10\B0)

    sfrac = 1.0  # fractional solar constant (def = 1.0)
    ktrans = 3.81  # transport coefficient (W m^(-2) K^(-1)) -  relaxation coeffient in Newtonian (linear) approximation of meridional heat transport divergence (def = 3.81)
    aice = 0.62  # albedo of snow and ice (def = 0.62)
    alw = 204  # constant part of linearized longwave radiation (W m^(-2)) (def = 204)
    blw = 2.17  # linear part of linearized longwave radiation (W m^(-2) K^(-1)) (def = 2.17)
    t0 = 15.0  # initial temperature (deg C) (def = 15)

    # Initialize time step
    dt = 30          # time step in days
    dt = dt * mpf(86400)    # time step in seconds
    year0 = 1.0        # initial model year
    model_time = year0  # current date time
    runlen = 100.0     # run length in years
    ndaymonth = 30.0  # number of days per month
    year_len = 12 * ndaymonth  # number of days in year

    ittmax = round(runlen * year_len * mpf(86400) / dt)  # maximum number of integer time counter (maximum number of timesteps)

    latzone = df_init['Zone']
    s = df_init['s']
    asfc = df_init['asfc']

    # ======== Running the model with the defined input variables =========
    df_EBM, tmean, yice_s, yice_n, jice_s, jice_n = EnergyBalanceModel1D(ittmax, grid_step, asfc, year0, year_len, dt,
                                                                         cocean, s, latzone, scon0, tzero, tcrit, t0,
                                                                         sfrac, ktrans, aice, alw, blw, to_print=True)


    # ======== Find "Glaciated Earth Edge Point" and plot graph of temperature vs sfrac,    =========
    #          going from non-glaciated Earth to glaciated Earth

    sfrac_range = [0.5, 1.8]
    step = 0.01

    sfrac_ice, df_sfrac = find_sfrac_ice(sfrac_range, step, ittmax, grid_step, asfc, year0, year_len, dt, cocean, s, latzone,
                                         scon0, tzero, tcrit,  t0, ktrans, aice, alw, blw)
    print("Fractional solar constant value at Glaciated Earth Edge Point is:", sfrac_ice)

    # df_sfrac.to_csv('sfrac.csv')  # save data to csv

    title = 'Final global mean temperature vs fraction of solar constant'
    xlabel = 'Fraction of solar constant (sfrac)'
    ylabel = 'tmean (deg C)'

    plot_estimated_line(df_sfrac, x='sfrac', y='tmean', vline_place=sfrac_ice,
                        title=title, xlabel=xlabel, ylabel=ylabel, grid=True, linewidth=1, color='blue', legend=False)

    plt.show()


    # ======== Find "Unglaciated Earth Edge Point" and plot graph of temperature vs sfrac,      ========
    #          going from glaciated Earth to non-glaciated Earth

    sfrac_range = [1.8, 0.5]
    step = 0.01

    sfrac_ice, df_sfrac = find_sfrac_ice(sfrac_range, step, ittmax, grid_step, asfc, year0, year_len, dt, cocean, s, latzone,
                                         scon0, tzero, tcrit,  t0, ktrans, aice, alw, blw)

    print("Fractional solar constant value at Unglaciated Earth Edge Point is:", sfrac_ice)

    # df_sfrac.to_csv('sfrac_inverse.csv')  # save data to csv

    title = 'Final global mean temperature vs fraction of solar constant'
    xlabel = 'Fraction of solar constant (sfrac)'
    ylabel = 'tmean (deg C)'

    plot_estimated_line(df_sfrac, x='sfrac', y='tmean', vline_place=sfrac_ice,
                        title=title, xlabel=xlabel, ylabel=ylabel, grid=True, linewidth=1, color='blue', legend=False)

    plt.show()


    #======== Find sfrac values for Glaciated Earth Edge Point, for different B (blw) values ========

    sfrac_range = [0.5, 1.1]
    step = 0.01
    blw_range = [1.0, 2.5]
    blw_step = 0.1

    blw_sfrac_ice = pd.DataFrame(columns=['blw', 'sfrac_ice'])
    for j in np.arange(blw_range[0], blw_range[1], step=blw_step):
        sfrac_ice, df_sfrac = find_sfrac_ice(sfrac_range, step, ittmax, grid_step, asfc, year0, year_len, dt, cocean, s,
                                             latzone, scon0, tzero, tcrit,  t0, ktrans, aice, alw, blw=j, run_all=False)
        blw_sfrac_ice = blw_sfrac_ice.append({'blw': j, 'sfrac_ice': sfrac_ice}, ignore_index=True)
    blw_sfrac_ice = blw_sfrac_ice.astype(float)
    blw_sfrac_ice.to_csv('blw_sfrac_ice.csv')

    title = 'Glaciated Earth Edge Point for different B values'
    xlabel, ylabel = 'B', 'sfrac for glaciated Earth'
    ax = plot_estimated_line(blw_sfrac_ice, x='blw', y='sfrac_ice',
                        title=title, xlabel=xlabel, ylabel=ylabel, grid=True, linewidth=1, color='blue')
    ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))
    plt.show()


    # ======== Climate for different values of Transport Coeff (ktrans) ========

    ktrans_range = [2.0, 5.0]
    ktrans_step = 0.1

    df_ktrans, df_ktrans_tmean = pd.DataFrame(), pd.DataFrame()

    for ktrans in np.arange(ktrans_range[0], ktrans_range[1], ktrans_step):
        df, tmean = EnergyBalanceModel1D(ittmax, grid_step, asfc, year0, year_len, dt, cocean,  s, latzone,
                                         scon0, tzero, tcrit, t0, sfrac, ktrans, aice, alw, blw, to_print=False)[0:2]

        df.insert(0, 'ktrans', ktrans)
        df_ktrans = df_ktrans.append(df, ignore_index=True)
        df_ktrans_tmean = df_ktrans_tmean.append({'ktrans': ktrans, 'tmean': tmean}, ignore_index=True)

    df_ktrans.loc[:, df.columns != 'Zone'] = df_ktrans.loc[:, df.columns != 'Zone'].astype(float)
    df_ktrans.to_csv('df_ktrans.csv')
    df_ktrans_tmean = df_ktrans_tmean.astype(float)
    df_ktrans_tmean.to_csv('df_ktrans_tmean.csv')

    title = 'Temperature for different values of transport coeff'
    xlabel, ylabel = 'ktrans', 'Temperature (deg C)'
    ax = plot_estimated_line(df_ktrans[df_ktrans.Zone == '80-90'], x='ktrans', y='Temp (deg C)',
                        title=title, xlabel=xlabel, ylabel=ylabel, grid=True, linewidth=1, color='blue')
    plot_estimated_line(df_ktrans[df_ktrans.Zone == ' 0-10'], x='ktrans', y='Temp (deg C)',
                             title=title, xlabel=xlabel, ylabel=ylabel, grid=True, linewidth=1, color='green', ax=ax)
    plot_estimated_line(df_ktrans_tmean, x='ktrans', y='tmean',
                             title=title, xlabel=xlabel, ylabel=ylabel, grid=True, linewidth=1, color='grey', ax=ax)

    # creating legend
    custom_lines = [plt.Line2D([0], [0], color='blue', lw=2),
                    plt.Line2D([0], [0], color='green', lw=2),
                    plt.Line2D([0], [0], color='grey', lw=2)]

    ax.legend(custom_lines, ['Poles', 'Equator', 'Mean'])

    plt.show()
    
