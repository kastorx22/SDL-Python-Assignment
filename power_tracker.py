#Kateryna
import random
#create empty list and starting index
nums=[]
i=0
flag=True
#while the condition is satisfied, random numbers are generated, power chosen, result is calculated and appended to the list
while flag:
    n=random.randint(1, 20)
    power=random.randint(2,3)
    num=n**power
    nums.append(num)
    print(f'loop{i+1} : {n}^{power}={num}')

    #done by Ariel, decision making to see if the loop needs to break
    if (nums[i] % nums[i - 1] == 0) and (i != 0) and (nums[i - 1] != 1):
        flag = False   
    else:
        i+=1

#Ariel
#sort numbers in the list for further processing
nums_sorted=sorted(nums)   
# output data
print(f"The largest result is {nums_sorted[-1]}")
print(f"The smallest result is {nums_sorted[0]}")
print(f"{nums[-1]} is divisible by {nums[-2]}.")
print(f"We completed {i+1} loops.")
