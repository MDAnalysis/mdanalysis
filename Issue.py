#issue easy inverse index array
from typing import List
import array
#method 01
List=input(list)
def inverse(list):
   
    res=list[::-1]
    print(res)
l=inverse(List)
print(l)

#method 02
def inverse_array(list):
    srr=[2,3,4,5,67,8]
    result=list(reversed(srr))
    print(result)
    
    
