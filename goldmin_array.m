function [fx,index]=goldmin_array(f)
f=-1*f;
f_c=f(1);
x_i=1;
for i=1:length(f)
    if f_c>f(i)
        f_c=f(i);
        x_i=i;
        i=i+1;
    end
end
index=x_i;
fx=-1*f_c;