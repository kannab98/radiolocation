function z = test_function()

file = dlmread('test.tsv');
x = file(:,1);
y = file(:,2);
z = x + 1i .* y;
end
