function [y,s]=iid(sigma,lambda,n)

ybar=lambda*sigma;
y=linspace(-ybar,ybar,n);
s=zeros(n,1);
for j=1:n
    s(j)=nchoosek(n-1,j-1);
end
s=s/2^(n-1);

end