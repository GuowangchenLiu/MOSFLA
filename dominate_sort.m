function x_frontvalue = dominate_sort(x,s,r)
functionvalue = x(:,s+1:s+r);
fnum=0;                                             %��ǰ�����ǰ������
cz=false(1,size(x,1));                  %��¼�����Ƿ��ѱ�������
frontvalue=zeros(size(cz));                         %ÿ�������ǰ������
[functionvalue_sorted,newsite]=sortrows(functionvalue); %���ݵ�һ�е���ֵ�������ƶ�ÿһ�У������һ�е���ֵ����ͬ�ģ��������ұȽ�
while ~all(cz)                                      %��ʼ�����ж�ÿ�������ǰ����,���øĽ���deductive sort
    fnum=fnum+1;
    d=cz;
    for i=1:size(functionvalue,1)
        if ~d(i)
            for j=i+1:size(functionvalue,1)
                if ~d(j)
                    k=1;
                    for m=2:size(functionvalue,2)
                        if functionvalue_sorted(i,m)>functionvalue_sorted(j,m)
                            k=0;
                            break
                        end
                    end
                    if k
                        d(j)=true;
                    end
                end
            end
            frontvalue(newsite(i))=fnum;
            cz(i)=true;
        end
    end
end
x_frontvalue=frontvalue';
 end
