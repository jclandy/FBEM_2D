function y = IF(condition, x1, x2)
% IF M-file for IF ternary operator
% Syntax: y = IF(condition, x1, x2)
%
% The elements of the output matrix will be set as:
% y(condition==true) = x1(condition==true)
% y(condition==false) = x2(condition==false)
%
% x1 and x2 can be function handles or matrices the same size as condition
%
% Example:
% x = -10:10;
% y = questioncolon(x==0, @(z) 1, @(z) sin(x(z))./x(z));

x1=x1.*ones(size(condition));
x2=x2.*ones(size(condition));

y = zeros(size(condition));
y(condition) = x1(condition);
y(~condition) = x2(~condition);

end