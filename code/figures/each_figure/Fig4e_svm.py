#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  2 23:25:55 2024

@author: kaitong
"""

from sklearn.svm import SVC
import pandas as pd

# Load the data
data = pd.read_csv('data_summary.csv')

# Fit the SVM model
svm_model = SVC(kernel='linear')
svm_model.fit(X=data[['Mean_volume', 'Mean_AR']], y=data['ColonyMorph'])

# Extract and calculate slope and intercept of the decision boundary
# Extract the coefficients
w = svm_model.coef_[0]
# The slope of the decision boundary
slope = -w[0] / w[1]
# Extract the intercept
b = svm_model.intercept_[0]
# The intercept of the decision boundary
intercept = -b / w[1]
# Print results
print(slope, intercept)


