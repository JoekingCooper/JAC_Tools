import numpy as np
class GridSearch:
    """
        This class is designed to propagate the errors on a function using a grid search method.
        
        It requires a function to be assigned to the class, as well as the input parameters and errors for that function.
    """
    def __init__(self,values_in,value_errors_in,function_in):
        """
            Initialize the class, assign the function to the class attributes and assign the variables to class parameters.
        """
        self.minimum_error=[]
        self.maximum_error=[]
        self.true_values=tuple(values_in)
        self.true_errors=value_errors_in
        self.temp_values=values_in
        self.function=function_in
    def errorpropagation(self,steps):
        """
            Decide which propagation method is needed depending on whether there are single or multiple variables.
        """
        if len(self.true_values)==1:
            self.true_values=self.true_values[0]
            return self.errorpropagation_single(steps)
        else:
            return self.errorpropagation_multiple(steps)
    def errorpropagation_multiple(self,steps):
        """
            Propagate the errors on a function for multiple inputs.
        """
        self.loop_rec(0,steps)
        return min(self.minimum_error),max(self.maximum_error),self.function(self.true_values)
    def loop_rec(self,n,steps):
        """
            The primary function.
            
            This recursive loop creates a grid of input parameters based on the number of dimensions and steps of the input.
        """
        if n <= len(self.true_errors)-1:
            ers=self.true_errors[n]
            x=np.linspace(self.true_values[n]+ers[0],self.true_values[n]+ers[1],steps)
            for i in x:
                self.temp_values[n]=i
                self.loop_rec( n + 1)
        else:
            R_temp=self.function(self.temp_values)
            R_true=self.function(self.true_values)
            if R_temp-R_true>0:
                self.maximum_error.append(R_temp-R_true)
            if R_temp-R_true<0:
                self.minimum_error.append(R_temp-R_true)
    def errorpropagation_single(self,steps):
        """
            Propagate the errors on a function for a single parameter.
        """
        ers=self.true_errors[0]
        x=np.linspace(np.array(self.true_values)+ers[0],np.array(self.true_values)+ers[1],steps)
        for i in x:
            temp=i
            R_temp=self.function(temp)
            R_true=self.function(self.true_values)
            if R_temp-R_true>0:
                self.maximum_error.append(R_temp-R_true)
            if R_temp-R_true<0:
                self.minimum_error.append(R_temp-R_true)
        R=self.function(self.true_values)
        return min(self.minimum_error),max(self.maximum_error),R
        
