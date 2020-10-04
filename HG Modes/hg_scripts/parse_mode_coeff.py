# Parse .txt's from Mathematica into Python lists

# 1. Mathematica solution 
mat_sol_arr = [None]*crossterm_order
for f in range(crossterm_order):
    file = open('1shift'+str(f)+".txt")

    line = file.read().replace("\n", " ")
    file.close()

#     print(line)
    print(f)
    mat_sol_arr[f] = line

mat_sol=mat_sol_arr[working]

# 2. Parse for Python

# parse solution from Mathematica FortranForm
# pars_sol=mc.mathematica(mat_sol)
# print(pars_sol)

str_pars_sol = str(mat_sol) #entire solution as string



# print(str_pars_sol.split(' '))
sol_list = str_pars_sol.split(" ") # split by whitespace



#remove '+', prepend '-'
for i in range(len(sol_list)):
    j = sol_list[i]
    if(j == '-'):
        sol_list[i+1]= '-' + sol_list[i+1]



# whitespace in Factorial(...) interferes w/ strip
## Example:
## (-8*a*b**3*p**4*z**3*Sqrt(Factorial(-4 + n)))/(3.*w**4)
## -> ["(z**3*Sqrt(Factorial(-4", "n)))/(3.*w**4)" ]
## 1) "(z**3*Sqrt(Factorial(-4" -> "(z**3*Sqrt(Factorial("
## 2) "n)))/(3.*w**4)" -> ")))/(3.*w**4)" 
## 3) -> ["(z**3*Sqrt(Factorial(n-4)))/(3.*w**4)" ]

#build list of terms, where operators +,-, whitespace deleted
sol_list[:] = [i for i in sol_list if i not in ['-', '+', '']] 

def parse_factorial(sol_list):
    for i in sol_list:
        if 'Factorial' in i:
            current_index = sol_list.index(i) #index with factorial
            next_index = current_index+1
            n_minus = 'n-' + i[-1]
            sol_list[current_index] = sol_list[current_index][:-2] + n_minus +  sol_list[next_index][1:]
            del sol_list[next_index]

parse_factorial(sol_list)

            
# print("Terms list:", sol_list)

#3. Separate terms list for x-coordinate dependence

#get highest n for x**n
rows, cols = (10, 10) 
sols_matrix = [ [ '' for i in range(rows) ] for j in range(cols) ]

#remove x**n and put term into sol_list[n]
for i in sol_list:
    x_order=0; p_order = 0
    #if x -> arr[n]
    #nonlinear
    if('x**' in i):
        exp_ind=i.find('x')+3 #get exp index
        x_order = int(i[exp_ind]) #set x order
    #no x -> arr[0]
    elif(not 'x' in i):
        x_order = 0
        #linear 
    else:
        x_order = 1
        
    #if x -> arr[n]
    #nonlinear
    if('p**' in i):
        exp_ind=i.find('p')+3 #get exp index
        p_order = int(i[exp_ind]) #set x order
    #no x -> arr[0]
    elif(not 'p' in i):
        p_order = 0
        #linear 
    else:
        p_order = 1
        
    sols_matrix[x_order][p_order] += '+' + i  
#     print(i,x_order,p_order)

    

    
# for row in range(rows):
#     for col in range(cols):
#         print(sols_matrix[row][col],row,col)

# print(sols_matrix)