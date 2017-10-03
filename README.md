# SVM_Framework
Support vector machines flexible framework  
We solve the unconstrained primal SVM formulation using representer theorem  
That is, we use different types of Kernels: linear, gaussian and Polynomial.  
Customized kernels can be easily added as well.  
The framework supports several features with several updates ongoing:  
1. Binary classification:  
-huberized loss minimization using newton optimization  
-Implementation of SGD using pegasos algorithm  
2. Multiclass classification:  
SVM hinge loss minimization using gradient descent  
Convergence time for this mode is slow. I'll provide faster implementations asap.  
3. Regression  
