
% ARO2021L Section 5
% Week 5 | Help/Hints
% Rohan Patel


% Extra. clear your Command Window with:
clc


% Extra. use ";" at the end of your code line to supress the output
'this displays'
'this does not';












% =========================================================================
% 1. Accessing MATLAB's documentation
%       "help" followed by the name of the MATLAB function
%       The MATLAB docs are VERY helpful and have examples!
help diff
help plot
help rand










% =========================================================================
% 2. Array Basics
%       Arrays can be 1-D : [1, 2, 5, -2, 1000]
%       Arrays can be 2-D : [1, 2; 5, -2; -222 1000]
row1Darray = rand(1,3);     % <--- row vector (or array) of height=1, length=3
col1Darray = rand(3,1);     % <--- column vector (or array) of height=3, length=1
example2Darray = rand(4);   % <--- 2-D array of height=4, length=4
exampleArray = 0.0:0.5:10;  % <-- 1-D array of values starting at 0 and going to 10 in increments of 0.5

% Array Indicies:
%   Remember: arrayName(row, col) <-- MATLAB starts from 1 not 0!
% lets say I want the first element in array "row1Darray"
firstElement = row1Darray(1);

% lets say I want the third element in the second row in "example2Darray"
ele2darray = example2Darray(2, 3);

% lets say I want the last element of the last row in "example2Darray"
lastlast = example2Darray(end, end);

% lets say I want the last 2 elements in row1Darray
lastThree = row1Darray(end-2 : end)     % <-- be careful with commas v. colons
% colons = selecting a range of values
% commas = sep. between rows and cols










% =========================================================================
% 2. Using diff() 
data1 = [1,2,3,4,5];
diff1 = diff(data1);
% If we are taking the difference between each element, the size of 
% the differenced array will never equal the size of the original array!












% =========================================================================
% 3. Plotting

xData = rand(100,1);
yData = rand(100,1);

plot(xData, yData)  % <--- plots array xData versus array yData

% Let's add axis labels and tell MATLAB to make the fontsize big
xlabel('My X-Data Name', 'FontSize', 20);
ylabel('My Y-Data Name', 'FontSize', 20);

% Add a title
title('this is my title')

% Add a x/y grid
grid on;

% Make my axis lengths equal
axis equal;












% =========================================================================
% 4. Functions

% Lets say I want to make a function that adds two numbers and subtracts
% them:


% Calling your function:
cat = 100;
dog = 329;
[catdog, somenamehere] = myFunction(cat,dog)
% cat = first input being passed into the function
% dog = second input being passed into the function
% catdog = first output (which corresponds to addedValues)
% somenamehere = second output (which corresponds to subtractedValues)


% Call my function again, and again!:
b = 029384;
u = 9024820938320;
[output1, output2] = myFunction(b,u);
[output11, output22] = myFunction(u,b);
[output111, output222] = myFunction(u,u);




% Defining your function:
function [addedValues, subtractedValues] = myFunction(a,b)
% a = first input
% b = second input
% myFunction = name of the function (which we use to call it in our code)
% addedValues = first output
% subtractedValues = second output


moreMathIcandohere = a-b;

addedValues = a+b;
subtractedValues = moreMathIcandohere;

end



