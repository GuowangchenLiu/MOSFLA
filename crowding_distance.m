function P = crowding_distance(X,s,r)
P=zeros(size(X,1),1);
front_number=max(X(:,s+r+1));
functionvalue=X(:,s+1:s+r);
for m=1:front_number
    popu=find(X(:,s+r+1)==m)';                                          %popu��¼��fnum+1�����ϵĸ�����
    distancevalue=zeros(size(popu));                                        %popu�������ӵ������
    fmax=max(functionvalue(popu,:),[],1);                                   %popuÿά�ϵ����ֵ
    fmin=min(functionvalue(popu,:),[],1);                                   %popuÿά�ϵ���Сֵ
    for i=1:size(functionvalue,2)                                           %��Ŀ�����ÿ��Ŀ����popu�������ӵ������
        [~,newsite]=sortrows(functionvalue(popu,i));
        distancevalue(newsite(1))=inf;
        distancevalue(newsite(end))=inf;
        for j=2:length(popu)-1
            distancevalue(newsite(j))=distancevalue(newsite(j))+(functionvalue(popu(newsite(j+1)),i)-functionvalue(popu(newsite(j-1)),i))/(fmax(i)-fmin(i));
        end
    end
    P(popu,1)=distancevalue';
end
end
