function f = evaluate_objective(x)
% Function to evaluate the objective functions is ZDT1
% x. x has the decision variables
f=zeros(size(x,1),2);           %�ϲ�����Ⱥ�ĸ�Ŀ�꺯��ֵ,�����������ZDT1
f(:,1)=x(:,1);                  %�����һάĿ�꺯��ֵ
g=1+9*sum(x(:,2:30),2)./29;
f(:,2)=g.*(1-(x(:,1)./g).^0.5); %����ڶ�άĿ�꺯��ֵ
end
