Note:

1. Each specie configuration file is unique to the selected equation files. 
   Users cannot use specie configuration file generated for one set of equations for another set of equations.

2. Specie unit only accepts ppbv or molec_cm3

3. For any rate expression that involves calculation, there must be a space before and after a function call.
   There cannot be any space between arguments for a function.
   For numeric operands (including numeric variables), there could be any number of spaces between them and operators.
   The first character in expression cannot be a parenthese.

   Examples of correctly formatted rate expressions:
       a) 0.21* ARR2(3.30e-11,-55,TEMP) +0.78* ARR2(2.15e-11,-110,TEMP)
       b) 1.44e-13*(1+(C_M/4.2e+19))
       c) k37(TEMP,C_M,C_H2O)
       d) 7.2e-34 *0.78084* C_M*0.20946*C_M * (TEMP / 300)**(-2.6)

   Examples of incorrectly formatted rate expressions:
       a) (0.21*ARR2( 3.30e-11,-55.0,TEMP) + 0.78*ARR2(2.15e-11,-110,TEMP))
       b) k37( TEMP, C_M, C_H2O )

4. Each quation file must contains at least one line of comments above all equations
   The comments must contains at least three commas to ensure separation of columns functions correctly

5. hv and duplication skipping are disabled in eqns.py

6. M, O2, and hv are not presented in KPP model