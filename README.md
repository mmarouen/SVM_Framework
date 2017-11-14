# SVM_Framework
Support vector machines flexible framework  
We solve the unconstrained primal SVM formulation using representer theorem  
The framework is flexible  
#1. Classifiers/regressors:  
NB: multinomial classification under dev  
LS: regression classifier using penalized least squared loss  
Softmax: Softmax classifier using cross entropy loss  
SVM:svm classifier using quadratic hinge loss  
#2. Optimization methods:  
BGD: gradient descent (batch)  
NGD: Newton-Raphson optimization (batch)  
CGD: conjugate gradient descent (batch)    
SGD: stochastic gradient descent (under development)  
#3. Kernels:  
gaussian: gaussian kernel  
linear: linear kernel  
poly:polynomial kernel  
