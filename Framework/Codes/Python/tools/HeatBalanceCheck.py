# -*- coding: utf-8 -*-
"""
# Energy Balance Check

The goal is here to

* Quantify the impact of HR at household level
* Quantify the impact of HR at ARA level

The difference may be explained by lateral connections
"""
from __future__ import print_function
import pandas as pd
from BinaryOutput import readFromBinary

if __name__=="__main__":
    c_water = 4179.6 ## J/L*K
    refT = 10 ## °C

    ## Energy Balance Residential (with/without HR)
    workingDir = 'C:/Users/hadengbr/04_Programming/tmp/SWMM-HEAT'
    # workingDir = "Q:/Abteilungsprojekte/eng/SWWData/SWMM-HEAT/Framework_template/Residential/Households/stitched"
    nodesDf = pd.read_excel('Q:/Abteilungsprojekte/eng/SWWData/SWMM-HEAT/Framework_template/Data/Residential.xlsx')
    print(nodesDf["Node"])

    totalHeat_noHR_dry = 0.0
    totalHeat_HR_dry = 0.0
    totalHeat_noHR_wet = 0.0
    totalHeat_HR_wet = 0.0

    for node in nodesDf["Node"]:
    # for node in ['5h']:
        print("node {}".format(node))
        houses = int(round(nodesDf[nodesDf['Node'] == node]['HousesPerNode']))

        if houses == 0:
            continue

        for house in range(1, houses + 1, 1):
            print("reading house {}...".format(house))
            df_noHR = readFromBinary("{}/Reference_20190302-20190316/{}/FlowTemp_{}_{}.bin".format(workingDir, node, node, house), ['flow', 'temperature']).resample('1S').bfill()
            df_HR = readFromBinary("{}/PHR_1_100_20190302-20190316/{}/FlowTemp_{}_{}.bin".format(workingDir, node, node, house), ['flow', 'temperature']).resample('1S').bfill()

            # print(df_noHR)

            df_noHR_dry = df_noHR.loc["2019-03-04 00:00:00":"2019-03-09 00:00:00"]
            df_HR_dry = df_HR.loc["2019-03-04 00:00:00":"2019-03-09 00:00:00"]
            totalHeat_noHR_dry += (df_noHR_dry['flow'] * (df_noHR_dry['temperature'] - refT)).sum() * c_water / 3600000  ## times 5 because of time resolution 5 seconds
            totalHeat_HR_dry += (df_HR_dry['flow'] * (df_HR_dry['temperature'] - refT)).sum() * c_water / 3600000  ## times 5 because of time resolution 5 seconds

            df_noHR_wet = df_noHR.loc["2019-03-11 00:00:00":"2019-03-16 00:00:00"]
            df_HR_wet = df_HR.loc["2019-03-11 00:00:00":"2019-03-16 00:00:00"]
            totalHeat_noHR_wet += (df_noHR_wet['flow'] * (df_noHR_wet['temperature'] - refT)).sum() * c_water / 3600000  ## times 5 because of time resolution 5 seconds
            totalHeat_HR_wet += (df_HR_wet['flow'] * (df_HR_wet['temperature'] - refT)).sum() * c_water / 3600000  ## times 5 because of time resolution 5 seconds

    print("no heat recovery total heat flow - dry: {} kWh".format(totalHeat_noHR_dry))
    print("heat recovery total heat flow - dry: {} kWh".format(totalHeat_HR_dry))
    print("no heat recovery total heat flow - wet: {} kWh".format(totalHeat_noHR_wet))
    print("heat recovery total heat flow - wet: {} kWh".format(totalHeat_HR_wet))
