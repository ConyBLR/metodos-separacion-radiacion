# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 12:25:35 2024

@author: Cony
"""

import numpy as np
import pandas as pd

class Engerer2:
    def __init__(self):
        self.C = 4.2336e-2
        self.beta0 = -3.7912
        self.beta1 = 7.5479
        self.beta2 = -1.0036e-2
        self.beta3 = 3.1480e-3
        self.beta4 = -5.3146
        self.beta5 = 1.7073
        
    def dhi(self, data):
        df = data[data['sza'] < 83]
        kt = np.where(df['eth'] == 0, np.nan, df['ghi'] / df['eth'])
        ktc = np.where(df['eth'] == 0, np.nan, df['ghics'] / df['eth'])
        delta_ktc = ktc - kt
        ghi = df['ghi']
        GHIcc = df['ghics']
        AST = df['hs']
        Z = df['sza']

        def kde(ghi, GHIcc):
            result = np.where(ghi == 0, np.nan, np.maximum(0, 1 - (GHIcc / ghi)))
            return result
        
        kde_values = kde(ghi, GHIcc)
        kde_filtered = np.where(kde_values <= 1, kde_values, np.nan) 
        
        df = df.replace([np.inf, -np.inf], np.nan)
        #kd_engerer2 = self.C + ((1 - self.C) / (1 +  np.exp(self.beta0 + self.beta2 * AST + self.beta3 * Z)))+ (self.beta5 * kde_filtered)
        #kd_engerer2 = self.C + ((1 - self.C) / (1 + np.exp(self.beta0 + self.beta1 * kt + self.beta2 * AST + self.beta3 * Z + self.beta4 * delta_ktc))) + (self.beta5 * kde_filtered)
        
        '''kd_engerer2 = np.where(np.isnan(kt) | np.isnan(delta_ktc), np.nan, self.C + ((1 - self.C) / (1 + np.exp(self.beta0 + self.beta1 * kt  + self.beta2 * AST + self.beta3 * Z + self.beta4 ))) + (self.beta5 * kde_filtered))
        
        
        df['kd_engerer2'] = kd_engerer2
        
        dhi_engerer2 = np.where((df['TOA'] == 0) | (df['GHImed'] == 0), np.nan, df['kd_engerer2'] * df['GHImed'])
        
        df['dhi_engerer2'] = dhi_engerer2'''
        en2 = ( self.beta0 + self.beta1 * kt + self.beta2 * AST + self.beta3 * Z + self.beta4 * delta_ktc)
        
        exponencial = []
        for i in en2:
            #e = np.exp(i)
            exponencial.append(i)
            
        kd_engerer2 =  self.C + ((1 - self.C) / (1 + np.exp(exponencial ))) + (self.beta5 * kde_filtered)
        
        df['kt'] = kt
        df['kde'] = kde_filtered
        df['ktc'] = ktc
        df['delta_ktc'] = delta_ktc
        df['kd_engerer2'] = kd_engerer2
        df['en2'] = en2
        
        df['dhiE2'] = df['ghi']*kd_engerer2

        return df


class Engerer3:
    def __init__(self):
        self.C = 0.1090
        self.beta0 = -2.0506e-2
        self.beta1 = 8.1249
        self.beta2 = -3.6234e-2
        self.beta3 = -4.1397e-2
        self.beta4 = -5.1045
        self.beta5 = 0
        
    def dhi(self, data):
        df = data[data['sza'] < 83]
        kt = np.where(df['eth'] == 0, np.nan, df['ghi'] / df['eth'])
        ktc = np.where(df['eth'] == 0, np.nan, df['ghics'] / df['eth'])
        delta_ktc = ktc - kt
        ghi = df['ghi']
        GHIcc = df['ghics']
        AST = df['hs']
        Z = df['sza']

        def kde(ghi, GHIcc):
            result = np.where(ghi == 0, np.nan, np.maximum(0, 1 - (GHIcc / ghi)))
            return result
        
        kde_values = kde(ghi, GHIcc)
        kde_filtered = np.where(kde_values <= 1, kde_values, np.nan) 
        
        df = df.replace([np.inf, -np.inf], np.nan)

        en3 = ( self.beta0 + self.beta1 * kt + self.beta2 * AST + self.beta3 * Z + self.beta4 * delta_ktc)
        
        exponencial = []
        for i in en3:
            #e = np.exp(i)
            exponencial.append(i)
            
        kd_engerer3 =  self.C + ((1 - self.C) / (1 + np.exp(exponencial ))) + (self.beta5 * kde_filtered)
        
        df['kt'] = kt
        df['kde'] = kde_filtered
        df['ktc'] = ktc
        df['delta_ktc'] = delta_ktc
        df['kd_engerer3'] = kd_engerer3
        df['en3'] = en3
        
        df['dhiE3'] = df['ghi']*kd_engerer3

        return df

class Engerer4:
    def __init__(self):
        self.C = 1.0562e-1
        self.beta0 = -4.1332
        self.beta1 = 8.2578
        self.beta2 = 1.0087e-2
        self.beta3 = 8.8801e-4
        self.beta4 = -4.9302
        self.beta5 = 4.4378e-1
        
    def dhi(self, data):
        df = data[data['sza'] < 83]
        kt = np.where(df['eth'] == 0, np.nan, df['ghi'] / df['eth'])
        ktc = np.where(df['eth'] == 0, np.nan, df['ghics'] / df['eth'])
        delta_ktc = ktc - kt
        ghi = df['ghi']
        GHIcc = df['ghics']
        AST = df['hs']
        Z = df['sza']

        def kde(ghi, GHIcc):
            result = np.where(ghi == 0, np.nan, np.maximum(0, 1 - (GHIcc / ghi)))
            return result
        
        kde_values = kde(ghi, GHIcc)
        kde_filtered = np.where(kde_values <= 1, kde_values, np.nan) 
        
        df = df.replace([np.inf, -np.inf], np.nan)
        en4 = ( self.beta0 + self.beta1 * kt + self.beta2 * AST + self.beta3 * Z + self.beta4 * delta_ktc)
        
        exponencial = []
        for i in en4:
            #e = np.exp(i)
            exponencial.append(i)
            
        kd_engerer4 =  self.C + ((1 - self.C) / (1 + np.exp(exponencial ))) + (self.beta5 * kde_filtered)
        
        df['kt'] = kt
        df['kde'] = kde_filtered
        df['ktc'] = ktc
        df['delta_ktc'] = delta_ktc
        df['kd_engerer4'] = kd_engerer4
        df['en4'] = en4
        
        df['dhiE4'] = df['ghi']*kd_engerer4

        return df

class Yang4:
    def __init__(self):
        #coeficientes de Yang2
        self.C = 0.0361
        self.beta0 = -0.5744
        self.beta1 = 4.3184
        self.beta2 = -0.0011
        self.beta3 = 0.0004
        self.beta4 = -4.7952
        self.beta5 = 1.4414
        self.beta6 = -2.8396
        
    def dhi(self, data, frecuencia = 60 ):
        df = data[data['sza'] < 83]
        kt = np.where(df['eth'] <= 0.1, np.nan, df['ghi'] / df['eth'])
        ktc = np.where(df['eth'] <= 0.1, np.nan, df['ghics'] / df['eth'])
        delta_ktc = ktc - kt
        ghi = df['ghi']
        GHIcc = df['ghics']
        AST = df['hs']
        Z = df['sza']
        # el kd e2 horario correponde a data['kd_e2']
        if frecuencia == 60: 
            kds60 =  data['kd_e2'] #data['kd_engerer2']#
        else:
            kds60 =  data['kd_engerer2']#data['kd_e2'] #


        def kde(ghi, GHIcc):
            result = np.where(ghi == 0, np.nan, np.maximum(0, 1 - (GHIcc / ghi)))
            return result
        
        kde_values = kde(ghi, GHIcc)
        kde_filtered = np.where(kde_values <= 1, kde_values, np.nan) 
        
        df = df.replace([np.inf, -np.inf], np.nan)
        
        
        argumento = (self.beta0 + self.beta1 * kt + self.beta2 * AST + self.beta3 * Z + self.beta4 * delta_ktc + self.beta6 * kds60)
        
        exponencial = []
        for i in argumento:
            #e = np.exp(i)
            exponencial.append(i)
        
        kd_yang4 = self.C + (1 - self.C) / (1 + np.exp(argumento)) + (self.beta5 * kde_filtered)
        
        df['kt'] = kt
        df['kde'] = kde_filtered
        df['ktc'] = ktc
        df['delta_ktc'] = delta_ktc
        df['kd_yang4'] = kd_yang4
        df['dhiY4'] = df['ghi']*kd_yang4
        
        return df


class Yang5: #-65.75,-24.75,0.0469174926392916,0.110920228662633,0.589199573592444,5
    def __init__(self):
        #clima 4
        self.C = 0.042971966#-0.010948960
        self.beta0 = -1.644372959#-0.921287782
        self.beta1 = 4.718078475#3.650149139
        self.beta2 = 0.014623824#0.007674206
        self.beta3 = 0.007453108#0.004936045
        self.beta4 = -3.352233222#-3.764652421
        self.beta5 = 1.251921688#1.364819177
        self.beta6 = -2.364771589#-2.118672147
        
    def dhi(self, data, frecuencia = 60 ):
        df = data[data['sza'] < 83]
        kt = np.where(df['eth'] <= 0.1, np.nan, df['ghi'] / df['eth'])
        ktc = np.where(df['eth'] <= 0.1, np.nan, df['ghics'] / df['eth'])
        delta_ktc = ktc - kt
        ghi = df['ghi']
        GHIcc = df['ghics']
        AST = df['hs']
        Z = df['sza']
        # el kd e2 horario correponde a data['kd_e2']
        if frecuencia == 60: 
            kds60 =  df['kd_e2'] #data['kd_engerer2']#
        else:
            kds60 =  df['kd_engerer2']#data['kd_e2'] #


        def kde(ghi, GHIcc):
            result = np.where(ghi == 0, np.nan, np.maximum(0, 1 - (GHIcc / ghi)))
            return result
        
        kde_values = kde(ghi, GHIcc)
        kde_filtered = np.where(kde_values <= 1, kde_values, np.nan) 
        
        df = df.replace([np.inf, -np.inf], np.nan)
        
        
        argumento = (self.beta0 + self.beta1 * kt + self.beta2 * AST + self.beta3 * Z + self.beta4 * delta_ktc + self.beta6 * kds60)
        
        exponencial = []
        for i in argumento:
            #e = np.exp(i)
            exponencial.append(i)
        
        kd_yang5 = self.C + (1 - self.C) / (1 + np.exp(argumento)) + (self.beta5 * kde_filtered)
        
        df['kt'] = kt
        df['kde'] = kde_filtered
        df['ktc'] = ktc
        df['delta_ktc'] = delta_ktc
        df['kd_yang5'] = kd_yang5
        df['dhiY5'] = df['ghi']*kd_yang5
        
        return df