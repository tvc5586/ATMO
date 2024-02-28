"""
Functions for calculating rate expression number using kinetic functions
"""

import math
import pandas as pd

# TODO: Replace eval with a fully functional compiler

__all__ = ("kinetic_function_update", "function_not_exist")

# ----------------------------- Reaction Functions -----------------------------

def ARR(A0, B0, C0, TEMP):
    """ARR kinetic function with respect to temperature"""
    return A0 * math.exp(-B0/TEMP ) * (TEMP/300)**C0

def ARR2(A0, B0, TEMP):
    """ARR2 kinetic function with respect to temperature"""
    return A0 * math.exp(-B0 / TEMP)

def TROE(k0_300K, n, kinf_300K, m, TEMP, C_M):
    """TROE kinetic function with respect to temperature and C_M"""
    zt_help = 300 / TEMP
    k0_T = k0_300K * zt_help**(n) * C_M
    kinf_T = kinf_300K * zt_help**(m)
    k_ratio = k0_T / kinf_T

    return k0_T / (1 + k_ratio) * 0.6**(1 / (1 + math.log10(k_ratio)**2))

def TROEE(A, B, k0_300K, n, kinf_300K, m, TEMP, C_M):
    """TROEE kinetic function with respect to temperature and C_M"""

    zt_help = 300 / TEMP
    k0_T = k0_300K * zt_help**(n) * C_M
    kinf_T = kinf_300K * zt_help**(m)
    k_ratio = k0_T / kinf_T
    troe = k0_T / (1 + k_ratio) * 0.6**(1 / (1 + math.log10(k_ratio)**2))

    return A * math.exp(-B / TEMP) * troe

def k46(TEMP, C_M):
    """ARR kinetic function with respect to temperature, and C_M"""
    k0 = 2.4e-14 * math.exp(460/TEMP)
    k2 = 2.7e-17 * math.exp(2199/TEMP)
    k3 = 6.5e-34 * math.exp(1335./TEMP) * C_M

    return k0 + k3 / (1 + k3 / k2)

def k37(TEMP, C_M, C_H2O):
    """k37 kinetic function with respect to temperature,C_M, and C_H2O"""
    KMT06 = 1+(1.40e-21 * math.exp(2200/TEMP)* C_H2O )
    k2 = 2.20e-13 * math.exp(600/TEMP)
    k3 = 1.90e-33* C_M * math.exp(980/TEMP)

    return KMT06 * ( k2 + k3 )

def THERMAL_T2(c,d,TEMP):
    return (TEMP ** 2) * c * math.exp(-d/TEMP)

# ----------------------------- Module Functions -----------------------------

def remove_closing_bracket(txt):
    """Reform the last object in a list, in this case the closing bracket"""

    Last_element = list(txt[-1])
    Last_element.pop(-1)
    txt[-1] = "".join(Last_element)

    return txt


def kinetic_functions(TEMP,C_M,C_H2O,table):
    """Give the value of the kinetic constant when provided with a kinetic constant and its parameters"""
    for i in range(len(table)):

        L = list(table[i])

        if L[0].isalpha():
            new_rate = 0
            
            txt = table[i].split("(", 1)
            txt = remove_closing_bracket(txt)
            
            if txt[0] == 'ARR2':
                parameters = txt[1].split(",")
                new_rate = ARR2(float(parameters[0]), float(parameters[1]), TEMP)
            elif txt[0] == 'ARR':
                parameters = txt[1].split(",")
                new_rate = ARR(float(parameters[0]), float(parameters[1]), float(parameters[2]), TEMP)
            elif txt[0] == 'TROEE':
                parameters = txt[1].split(",")
                new_rate = TROEE(float(parameters[0]), float(parameters[1]), float(parameters[2]), float(parameters[3]),
                                       float(parameters[4]), float(parameters[5]), TEMP, C_M)
            elif txt[0] == 'TROE':
                parameters = txt[1].split(",")
                new_rate = TROE(float(parameters[0]), float(parameters[1]), float(parameters[2]), float(parameters[3]), TEMP, C_M)
            elif txt[0] == 'k46':
                parameters = txt[1].split(",")
                new_rate = k46(TEMP, C_M)
            elif txt[0] == 'k37':
                parameters = txt[1].split(",")
                new_rate = k37(TEMP, C_M, C_H2O)
            elif txt[0] == 'THERMAL_T2':
                parameters = txt[1].split(",")
                new_rate = THERMAL_T2(float(parameters[0]), float(parameters[1]), TEMP)
            elif txt[0] == 'TEMP' or txt[0] == 'C_M' or txt[0] == 'C_H2O':
                continue
            
            table[i] = str(new_rate)

def kinetic_functions_check(table):
    """Check if functions in rate expressions are implemented"""
    for i in range(len(table)):

        L = list(table[i])

        if L[0].isalpha():
            txt = table[i].split("(", 1)
            txt = remove_closing_bracket(txt)
            
            if (txt[0] == 'ARR2' or 
                txt[0] == 'ARR' or 
                txt[0] == 'TROEE' or 
                txt[0] == 'TROE' or 
                txt[0] == 'k46' or
                txt[0] == 'k37' or
                txt[0] == 'THERMAL_T2' or
                txt[0] == 'TEMP' or 
                txt[0] == 'C_M' or 
                txt[0] == 'C_H2O' or
                txt[0] == 'math.exp'):
                continue
            else:
                raise NotImplementedError

def kinetic_function_update(list_of_filenames, TEMP, C_M, C_H2O, pid):
    # WARNING: Current implementation is UNSAFE due to the use of "eval" function

    """Run an update of the kinetic constants with temperature, C_M and C_H2O.

    Parameters
    ----------
    list_of_filenames: list
        List of paths for equation files
        
    TEMP: float
        Temperature (in K) of the environment
        
    C_M: float
        
    C_H2O: float
    
    pid: int
    	Current process id
    
    Returns
    -------
    list_updated_files: list
        List of paths for updated equation files
    """
    
    """Reads the background_eqn file into a panda dateframe"""
    list_updated_files = []
    List_of_defined_functions = ['ARR2', 'ARR','TROEE','TROE','k46','k37','THERMAL_T2']

    for path in list_of_filenames:
        file_path = path
        eqn_file = pd.read_csv(file_path)
        (rows, column) = eqn_file.shape

        """For each reaction, reads the rate constant expression and calculate it to give one value"""
        for row in range(rows):
            first_term = eqn_file.iloc[row,0]
            first_character = list(first_term)
            
            if first_character[0] != '#':   ###This skips the row if the first cell has a #
                rate = eqn_file.iloc[row,2]
                table = rate.split(" ")
                
                # If rate expression contains a function
                for func in List_of_defined_functions:
                  if func in ''.join(table):
                    kinetic_functions(TEMP,C_M,C_H2O,table)
                
                rate = ''.join(table)
                
                # Replace leftover C_Ms, C_H2Os, and TEMPs with actual value
                rate = rate.replace("C_M", str(C_M))
                rate = rate.replace("C_H2O", str(C_H2O))
                rate = rate.replace("TEMP", str(TEMP))
                
                # Convert rate expression string into a number using eval (NOT SAFT)
                eqn_file.iloc[row, 2] = eval(rate) # eval is NOT SAFE to use. The ideal way is to implement a full compiler

        """Copy the new updated rate constant into another excel file"""
        Updated_file_name = str(path[:-4]) + "_updated_" + str(pid) + ".csv"
        eqn_file.to_csv(Updated_file_name, index = False, header = True)
        list_updated_files.append(Updated_file_name)

    return list_updated_files

def function_not_exist(list_of_filenames):
    """Check if functions the user chooses to use is included in the program

    Parameters
    ----------
    list_of_filenames: list
        List of paths for equation files
    """
    
    for path in list_of_filenames:                    ### For each excel file
        file_path = path
        eqn_file = pd.read_csv(file_path)
        (rows, column) = eqn_file.shape
        
        for row in range(rows):                  ### For each row rate expression
            first_term = eqn_file.iloc[row,0]
            first_character = list(first_term)
            
            if first_character[0] != '#':       ###This skips the row if the first cell has a #
              rate = eqn_file.iloc[row, 2]
              
              if "C_M" in rate or "C_H2O" in rate or "TEMP" in rate: # Replace C_M, C_H2O, or TEMP in rate expression with numeric placeholder
                rate = rate.replace("C_M", "0")
                rate = rate.replace("C_H2O", "0")
                rate = rate.replace("TEMP", "0")
              
              table = rate.split(" ")  
              
              # Check if functions are implemented
              kinetic_functions_check(table)

# ----------------------------- Test Code -----------------------------

if __name__ == "__main__":
    list_of_filenames = ['racm.csv']
    TEMP = 50000
    C_M = 60000
    C_H2O = 33000
    print(function_not_exist(list_of_filenames))
    new_file_paths = kinetic_function_update(list_of_filenames,TEMP,C_M,C_H2O, 000)
    print(new_file_paths)
    
