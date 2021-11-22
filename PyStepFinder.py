# -*- coding: utf-8 -*-
"""
Created on Sun Nov 21 14:21:21 2021

@author: glham
"""

import numpy
import sys
from scipy.stats import f

class chisq_functions: #each function computes the UNWEIGHTED sum of squared residualsfrom that statistic
    ##works with multiple arrays, easier to generalize if the splitting is done outside the function
    def mean(arrays): #works with multiple arrays, easier to generalize if the splitting is done outside the function
        return sum([len(array)*(sum([x**2. for x in array])/len(array)-(sum(array)/len(array))**2.) for array in arrays])
    def mean_no_n(arrays): #works with multiple arrays, easier to generalize if the splitting is done outside the function
    # variance is (X-<X>)^2/n, so each observable in the array minus the mean of that array, the times n gets rid of the unneeded division
        return sum([len(array)*(sum([x**2. for x in array])/len(array)-(sum(array)/len(array))**2.) for array in arrays])
    def mean_w_weight(arrays): #works with multiple arrays, easier to generalize if the splitting is done outside the function
        return sum([len(array)*(sum([x**2. for x in array])/len(array)-(sum(array)/len(array))**2.) for array in arrays])/(numpy.sum([len(array) for array in arrays])-len(arrays)) #weighted chisq statistic
    def standard_deviation(arrays):
        return sum([len(array)*numpy.log(sum([x**2. for x in array])/len(array)-(sum(array)/len(array))**2.) for array in arrays])
    def mean_reduced(arrays): #reduced chisq. number of means to estimate determines
        return chisq_functions.mean(arrays)/len(arrays)
    def mean_no_n_reduced(arrays):
        return chisq_functions.mean_no_n(arrays)/len(arrays)
    def mean_w_weight_reduced(arrays):
        return chisq_functions.mean_w_weight(arrays)/len(arrays)
    
    string={}
    string['mean']=mean
    string['mean_no_n']=mean_no_n
    string['mean_w_weight']=mean_w_weight
    string['mean_reduced']=mean_reduced
    string['mean_no_n_reduced']=mean_no_n_reduced
    string['mean_w_weight_reduced']=mean_w_weight_reduced
    string['standard_deviation']=standard_deviation

def BIC(n,k,RSS): #bayesian information criterion in terms of the RSS
    #n=number of data points
    #L=RSS for statistic
    #k=number of parameters
    return n*numpy.log(RSS/n)+k*numpy.log(n)

    
def findchangepointsopt(y, x=False, statistic='mean', min_length=1, min_separation=0, window=0, threshold=1.,threshold_mode='difference',MaxChanges=numpy.Inf,outputstyle='indices', optimize=None):
    #choose point, divide signal, computes estimate of signal parameter for each section, measures deviation for each point
    #adds deviations section by section to find error, minimizes this by varying division point
    #threshold is the threshold of chi squared change required to count addition of new cut. each cut, individually, must meat this criteria
        #once threshold unsatisfied, adjust locations of current cuts to optimize if optimization==True
        #retry addition of new cuts after optimization. if succeed, do so and continue until threshold not met, readjust, etc.
    #if threshold unmet on first pass, exit and return cuts
    #approach is pretty brute force so pretty slow. ¯\_(ツ)_/¯
    statistics=['mean','rms','variance','slope','mean_w_weight','mean_no_n']
    if statistic not in statistics:
        sys.exit('You ave selected an invalid parameter. This command only supports the following: '+ str(statistics))
    y=[y] #embedding initial list into a list of lists for easy storage of sub-list representing segments
    y_split=y #initialize the list that stores the output list of lists
    fin=False #flag to check if the whole run is complete based on chisq criteria
    add=True
    while not fin:
        itchisq=chisq_functions.string['{0}'.format(statistic)](y_split) 
        #chisq at current level of iteration through arrays. 
        #calculated as mean of statistic against mean of whole subarray. Sums chisq across all subarrays
        while add:
            additions=[] #indices of points to split
            for i in range(len(y_split)):
                #looks through arrays a level at a time: check all current arrays, make splits, check all new arrays, etc.
                candidate=[numpy.Inf,i,1] #array of chisq, subarray, and point of candidate split
                segchisq=chisq_functions.string['{0}'.format(statistic)]([y_split[i]])  #checks chisq of current subarray
                if len(y_split[i])<=2*min_length: #makes sure the split arrays adhere to minimum length requirement before searching to save effort
                    continue
                for k in range(min_length,len(y_split[i])-min_length): #only search over range that preserve the minimum length of a segment
                    candchisq=chisq_functions.string['{0}'.format(statistic)]([y_split[i][:k]]+[y_split[i][k:]]) #check the chisq of segments following split at point k
                    if candchisq<candidate[0] and segchisq!=0.: #if this candidate beats the curret top contender, chooe it
                    #ignore segchisq=0 b/c that indicates the subarray is single-valued in terms of the test statistic
                    #ie it cannot be split further as a single test statistic value adequately explains the dataset 
                        candidate=[candchisq,i,k]
                if candidate[0]==numpy.Inf:
                    continue
                if threshold_mode=='difference':
                    #check if adding a split improves the chisq by more than a given difference
                    if segchisq-candidate[0] > threshold:
                        additions.append(candidate[1:]) #add a list of the subarray index and split point in that array
                elif threshold_mode=='ratio': #check if ratiometric chisq improvement meet threshold
                    if segchisq-candidate[0] > threshold*segchisq or segchisq==numpy.Inf:
                        additions.append(candidate[1:])
                elif threshold_mode=='f-test': #check f-test
                    if (f.cdf(candidate[0]/segchisq,2,1) < threshold or segchisq==numpy.Inf):
                        additions.append(candidate[1:])
                elif threshold_mode=='BIC':
                    if BIC(len(y_split[i]),1,segchisq)-BIC(len(y_split[i]),2,candidate[0]) > threshold:
                        additions.append(candidate[1:])
            if additions==[]:
                add=False
            else:
                #print(additions)
                temp=[]
                additions.append([numpy.Inf,numpy.Inf])
                additions=numpy.array(additions)
                for i in range(len(y_split)):
                    if i in additions[:,0]:
                        #print(additions[additions[:,0]==i,1][0])
                        temp.append(y_split[i][:int(additions[additions[:,0]==i,1][0])])
                        temp.append(y_split[i][int(additions[additions[:,0]==i,1][0]):])
                    else:
                        temp.append(y_split[i])
                y_split=temp
        currchisq=chisq_functions.string['{0}'.format(statistic)](y_split) #new global chisq
        if itchisq==currchisq:    #if no changes were made, exit everything
            fin=True
        elif itchisq>currchisq:   #if improvement was made, go back for another round
            add=True
        if itchisq<currchisq:   #if results became worse, uh-oh
            print("something went wrong and chisq increased")
            fin=True
    count=0
    indices=[]
    for i in y_split: #generate list of split indices in terms of original array
        count+=len(i)
        indices.append(count)
    #print(arrays)
    if outputstyle in 'indices':
        return indices[:-1]
    if outputstyle in 'complete': #summary statistics and useable outputs
        results={}
        results['changepoints']=indices[:-1]
        results['segmented arrays']=y_split
        results['segment means']=[numpy.mean(y_split) for y_split in y_split]
        results['segment lengths']=[len(y_split) for y_split in y_split]
        results['plottable means']=[]
        for i in range(len(results['segment means'])):
            results['plottable means'].extend(numpy.ones(len(y_split[i]))*results['segment means'][i])
        results['initial chi-squared']=chisq_functions.string['{0}'.format(statistic)](y)
        results['final chi-squared']=chisq_functions.string['{0}'.format(statistic)](y_split)
        return results                        
                    
                
                        
    
    
    
    
    
    
    
    
    
    
    
    
    
    