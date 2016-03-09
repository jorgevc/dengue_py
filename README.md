# dengue_py
------------
Dependecies 
----------
-numpy
-scipy

---------
Run the simulations
---------

python ./dyn_mosquitos.py
This will show the incidence and mosqito population in a plot.
The comparison with medium, high and low temperature will be shown.

python ./correlation.py
The correlation between incidence and mosquito population with different time lags will also be generated.
The maximal correlation and time lag will be displayed

The main functions to build a custom symulation are in the library "find_best_treshold_lib.py" 
described below
------------------------------
Important Files and functions
------------------------------
-datos_plot_lib.py:
Library that reads the data from disk and make interpotations for missing data.

-mosquitos_lib.py:
Library with utility functions of the model. The parameters of the model are stored here

-find_best_treshold_lib.py:
Library with utility functions to test different values of parameters in recurrent way.
Main functions:
    mosquitos(region,data): -> weekly_population
      region: int index of the region
      data: data object
      weekly_mosquitos: array of weakly mosquito population 
      Return the mosquito population with varing temperature 
    mosquitos_diapause(region,data):
      region: int index of the region
      data: data object
      weekly_mosquitos: array of weakly mosquito population 
      Return the mosquito population with varing temperature, precipitation and diapause
    mosquitos_precipitation(region,data):
      region: int index of the region
      data: data object
      weekly_mosquitos: array of weakly mosquito population 
      Return the mosquito population with varing temperature and precipitation
    mosquitos_precipitation_sinTemp(region,data):
      region: int index of the region
      data: data object
      weekly_mosquitos: array of weakly mosquito population 
      Return the mosquito population with varing temperature

---------
Licence
---------
/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 */
