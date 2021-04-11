%function NSGAII()
clear all;close all;clc;format compact;tic;hold on
    
%---��ʼ��/�����趨
    generations=100;                                %��������
    popnum=100;                                     %��Ⱥ��С(��Ϊż��)
    poplength=30;                                   %���峤��
    % tepmat ����Ϊ�ѵ������� B=repmat(A,m,n)
    minvalue=repmat(zeros(1,poplength),popnum,1);   %������Сֵ
    maxvalue=repmat(ones(1,poplength),popnum,1);    %�������ֵ    
    population=rand(popnum,poplength).*(maxvalue-minvalue)+minvalue;    %�����µĳ�ʼ��Ⱥ
    
%---��ʼ��������
    for gene=1:generations                      %��ʼ����
        
%-------���� 
        newpopulation=zeros(popnum,poplength);  %�Ӵ���Ⱥ
        for i=1:popnum/2                        %��������Ӵ�
            k=randperm(popnum);                 %����Ⱥ�����ѡ��������ĸ,�����ö�������������
            beta=(-1).^round(rand(1,poplength)).*abs(randn(1,poplength))*1.481; %������̬�ֲ�������������Ӵ�
            newpopulation(i*2-1,:)=(population(k(1),:)+population(k(2),:))/2+beta.*(population(k(1),:)-population(k(2),:))./2;    %������һ���Ӵ�           
            newpopulation(i*2,:)=(population(k(1),:)+population(k(2),:))/2-beta.*(population(k(1),:)-population(k(2),:))./2;      %�����ڶ����Ӵ�
        end
        
%-------����        
        k=rand(size(newpopulation));    %���ѡ��Ҫ����Ļ���λ
        miu=rand(size(newpopulation));  %���ö���ʽ����
        temp=k<1/poplength & miu<0.5;   %Ҫ����Ļ���λ
        newpopulation(temp)=newpopulation(temp)+(maxvalue(temp)-minvalue(temp)).*((2.*miu(temp)+(1-2.*miu(temp)).*(1-(newpopulation(temp)-minvalue(temp))./(maxvalue(temp)-minvalue(temp))).^21).^(1/21)-1);        %�������һ
        newpopulation(temp)=newpopulation(temp)+(maxvalue(temp)-minvalue(temp)).*(1-(2.*(1-miu(temp))+2.*(miu(temp)-0.5).*(1-(maxvalue(temp)-newpopulation(temp))./(maxvalue(temp)-minvalue(temp))).^21).^(1/21));  %���������
        
%-------Խ�紦��/��Ⱥ�ϲ�        
        newpopulation(newpopulation>maxvalue)=maxvalue(newpopulation>maxvalue); %�Ӵ�Խ�Ͻ紦��
        newpopulation(newpopulation<minvalue)=minvalue(newpopulation<minvalue); %�Ӵ�Խ�½紦��
        newpopulation=[population;newpopulation];                               %�ϲ�������Ⱥ
        
%-------����Ŀ�꺯��ֵ        
        functionvalue=zeros(size(newpopulation,1),2);           %�ϲ�����Ⱥ�ĸ�Ŀ�꺯��ֵ,�����������ZDT1
        functionvalue(:,1)=newpopulation(:,1);                  %�����һάĿ�꺯��ֵ
        g=1+9*sum(newpopulation(:,2:poplength),2)./(poplength-1);
        functionvalue(:,2)=g.*(1-(newpopulation(:,1)./g).^0.5); %����ڶ�άĿ�꺯��ֵ
        
%-------��֧������        
        fnum=0;                                             %��ǰ�����ǰ������
        cz=false(1,size(functionvalue,1));                  %��¼�����Ƿ��ѱ�������
        frontvalue=zeros(size(cz));                         %ÿ�������ǰ������
        [functionvalue_sorted,newsite]=sortrows(functionvalue);    %����Ⱥ����һάĿ��ֵ��С��������
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
        
%-------����ӵ������/ѡ����һ������        
        fnum=0;                                                                 %��ǰǰ����
        while numel(frontvalue,frontvalue<=fnum+1)<=popnum                      %�ж�ǰ���ٸ���ĸ�������ȫ������һ����Ⱥ
            fnum=fnum+1;
        end        
        newnum=numel(frontvalue,frontvalue<=fnum);                              %ǰfnum����ĸ�����
        population(1:newnum,:)=newpopulation(frontvalue<=fnum,:);               %��ǰfnum����ĸ��帴������һ��                       
        popu=find(frontvalue==fnum+1);                                          %popu��¼��fnum+1�����ϵĸ�����
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
        popu=-sortrows(-[distancevalue;popu]')';                                %��ӵ�����뽵�������fnum+1�����ϵĸ���
        population(newnum+1:popnum,:)=newpopulation(popu(2,1:popnum-newnum),:);	%����fnum+1������ӵ������ϴ��ǰpopnum-newnum�����帴������һ��        
    end

%---�������    
    fprintf('�����,��ʱ%4s��\n',num2str(toc));          %�������պ�ʱ
    output=sortrows(functionvalue(frontvalue==1,:));    %���ս��:��Ⱥ�з�֧���ĺ���ֵ
    plot(output(:,1),output(:,2),'*b');                 %��ͼ
    axis([0,1,0,1]);xlabel('F_1');ylabel('F_2');title('ZDT1')
%end
