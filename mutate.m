function x_mutation=mutate(xg)
x_mutation=xg;
select=ones(1,4);
select=randperm(size(xg,2),3);
select=sort(select);
delta=select(2)-select(1);
if select(3)+delta>size(xg,2)
    select(4)=select(3)-delta;
    x_mutation(select(1):select(2))=x_mutation(select(4):select(3));
else
    select(4)=select(3)+delta;
     x_mutation(select(1):select(2))=x_mutation(select(3):select(4));
end
