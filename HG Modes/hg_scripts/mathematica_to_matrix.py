#######################
# Iteratively parse mathematica files of n shift m tilt into matrix of dim[n x m] and write to file
#######################
## These files all have a single string representing mode coupling coefficients from respective expansion orders
## Filename follows "n_shift_m_tilt.txt"
#######################

from ast import literal_eval

def parse_math_for_python(max_shift, max_tilt):
    sol_list = math_string(max_shift,max_tilt).split(" ") # split math_string into list by whitespace
    
    for i in range(len(sol_list)): # prepend + and - in list as signs of the following term
        j = sol_list[i]
        if(j == '-'):
            sol_list[i+1]= '-' + sol_list[i+1]
    
    sol_list[:] = [i for i in sol_list if i not in ['-', '+', '']] # remove the isolated signs and empties
    
    sol_list = parse_factorial(sol_list) # parses factorial for Python eval function and returns the solutions list
    
    x_p_matrix = parse_x_p(max_shift,max_tilt,sol_list)
    
    return(x_p_matrix)
    
#1. Open each file and get its string
def math_string(max_shift,max_tilt):

    mathematica_dir = "inputs_mathematica/"
        
    filename = str(max_shift)+"_shift"+str(max_tilt)+"_tilt.txt"
    file = open(mathematica_dir+filename)
    mat_sol = file.read().replace("\n", " ")
    file.close()
            
    return(mat_sol)    

###############
    ## whitespace in Factorial(...) interferes w/ strip
    ## Example:
    ## (-8*a*b**3*p**4*z**3*Sqrt(Factorial(-4 + n)))/(3.*w**4)
    ## -> ["(z**3*Sqrt(Factorial(-4", "n)))/(3.*w**4)" ]
    ## 1) "(z**3*Sqrt(Factorial(-4" -> "(z**3*Sqrt(Factorial("
    ## 2) "n)))/(3.*w**4)" -> ")))/(3.*w**4)" 
    ## 3) -> ["(z**3*Sqrt(Factorial(n-4)))/(3.*w**4)" ]
###############

# 2. Fix factorial terms so they aren't split from original function (not generalized to double digit exp. orders)
def parse_factorial(sol_list):
    for i in sol_list:
        if 'Factorial' in i:
            current_index = sol_list.index(i) #index with factorial
            next_index = current_index+1
            n_minus = 'n-' + i[-1]
            sol_list[current_index] = sol_list[current_index][:-2] + n_minus +  sol_list[next_index][1:]
            del sol_list[next_index]
            
    return(sol_list)


#3. Separate terms list into matrix by x-coordinate and p dependence
def parse_x_p(max_shift,max_tilt,sol_list):
    crossterm_order = max_shift+max_tilt
    rows, cols = (crossterm_order+1, crossterm_order+1) 
    x_p_matrix = [ [ '' for i in range(rows) ] for j in range(cols) ]

    #remove x**n and put term into sol_list[n]; also sort p_ordder
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

        x_p_matrix[x_order][p_order] += '+' + i  
        
    return(x_p_matrix)

#get dictionary from file of x p matrices
def dict_from_xpfile(filename):
    file = open("inputs/"+filename, "r")
    contents = file.read()
    dictionary = literal_eval(contents)
    file.close()
    return(dictionary)