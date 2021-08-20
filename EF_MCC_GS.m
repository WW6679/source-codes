function [aa]=EF_GS_SOCP(n,N_n)
% % 拉普拉斯拟合
It=50;
epsilon=1e-2;
ww = -ones(1,N_n);
w = -ww/sum(ww);
for it=1:It
% % ==============================================================求解参数向量a
if(it<3)
cvx_clear
cvx_begin
cvx_solver sedumi
cvx_precision best
cvx_quiet(1)
variable a(6,1);
variable gama(N_n,1);
minimize (sum(-w*gama))
subject to
    (n'*a).^2<=gama;
% for i=1:N_n
%   norm([2*n(:,i)'*a;gama(i)-1])<=gama(i)+1;
% end
    norm([a(2)-epsilon;a(1)-a(3)])<=a(1)+a(3);
    a(2)<=0;
cvx_end
target=sum(-w*gama);
% tar=sum(exp(-abs(a'*n)/sigma));
cvx_clear
cvx_begin 
cvx_solver sedumi
cvx_precision best
cvx_quiet(1)
variable a1(6,1);
variable gama(N_n,1);
minimize sum(-w*gama)
subject to
%     abs(a1'*n)'<=gama;
     (n'*a1).^2<=gama;
% for i=1:N_n
%   norm([2*n(:,i)'*a1;gama(i)-1])<=gama(i)+1;
% end
    norm([a1(2)+epsilon;a1(1)-a1(3)])<=a1(1)+a1(3);
    a1(2)>=0;
cvx_end
target1=sum(-w*gama);
storetarget(1,it)=target;
storetarget(2,it)=target1;
% tar1=sum(exp(-abs(a1'*n)/sigma));
% storetar(1,it)=tar;
% storetar(2,it)=tar1;
    if(target<=target1)
        Ma=a;
    else
        Ma=a1;
    end
if it ==1
    aa = Ma;
else
    aa0 = aa;
    aa = Ma;
    if norm(aa-aa0)<1e-5
        break;
    end
end
else
    if(aa(2)<0)
        cvx_clear
        cvx_begin 
        cvx_solver sedumi
        cvx_precision best
        cvx_quiet(1)
        variable a(6,1);
        variable gama(N_n,1);
        minimize sum(-w*gama)
        subject to
%             abs(a'*n)'<=gama;
         (n'*a).^2<=gama;
% for i=1:N_n
%   norm([2*n(:,i)'*a;gama(i)-1])<=gama(i)+1;
% end
            norm([a(2)-epsilon;a(1)-a(3)])<=a(1)+a(3);
            a(2)<=0;
        cvx_end
        target=sum(-w*gama);
        storetarget(1,it)=target;
%         tar=sum(exp(-abs(a'*n)/sigma));
%         storetar(1,it)=tar;
    else
        cvx_clear
        cvx_begin 
        cvx_solver sedumi
        cvx_precision best
        cvx_quiet(1)
        variable a(6,1);
        variable gama(N_n,1);
        minimize sum(-w*gama)
        subject to
%             abs(a'*n)'<=gama;
        (n'*a).^2<=gama;
% for i=1:N_n
%   norm([2*n(:,i)'*a;gama(i)-1])<=gama(i)+1;
% end
            norm([a(2)+epsilon;a(1)-a(3)])<=a(1)+a(3);
            a(2)>=0;
        cvx_end
        target1=sum(-w*gama);
        storetarget(2,it)=target1;
%         tar=sum(exp(-abs(a'*n)/sigma));
%         storetar(2,it)=tar1;
    end
    aa0 = aa;
    aa = a;
    if norm(aa-aa0)<1e-5
        break;
    end
end
aaa2(:,it)=aa;
% % ==============================================================更新权重 
data_points = aa.'*n;
sigma_E = std(data_points);
X = sort(data_points,'ascend');
Q1=X(floor((N_n)/4));
Q3=X(floor(3*(N_n)/4));
r = (Q3-Q1)/1.34;
sigma = 1.06*min(sigma_E,r)*N_n^(-0.2);
ww = -exp(-(data_points).^2/(2*sigma^2));
w = -ww/sum(ww);
end
aaa2(:,it)=aa;
itLPit=it;
aa_wokmeans=aa;
% % data_points = aa.'*n;
% % w1=log(1+log(1+log(1+abs(data_points))));
% % IDX=kmeans(w1',2);
% % % %%====================================
% % me1=sum(abs(w(IDX==1))/sum(IDX==1));
% % me2=sum(abs(w(IDX==2))/sum(IDX==2));
% % if(me2<me1)
% % %     disp('2表示野值');
% %     outlier_number=sum(IDX==2);
% %     for i=1:N_n
% %         if(IDX(i)==2)
% %             plot(xx(i),yy(i),'b+');hold on;
% %         end
% %     end
% % else
% % %     disp('1表示野值');
% %     outlier_number=sum(IDX==1);
% % for i=1:N_n
% %     if(IDX(i)==1)
% %         plot(xx(i),yy(i),'b+');hold on;
% %     end
% % end
% % end
% % if((me2>me1&&sum(IDX==1)<=60)||(me2<me1&&sum(IDX==2)<=60))
% % %     disp('成功');
% %     O=1;
% % else
% % %     disp('失败');
% %     O=0;
% % end
% % if(O==0)
% %     for ki=1:10
% %     ww = -abs(rand(1,N_n));
% %     w = -ww/sum(ww);
% %     for it=1:It
% %         % % ==============================================================求解参数向量a
% %         if(it<3)
% %         cvx_clear
% %         cvx_begin
% %         cvx_solver sedumi
% %         cvx_precision best
% %         cvx_quiet(1)
% %         variable a(6,1);
% %         variable gama(N_n,1);
% %         minimize sum(-w*gama)
% %         subject to
% % %             abs(a'*n)'<=gama;
% %        (n'*a).^2<=gama;
% %             norm([a(2)-epsilon;a(1)-a(3)])<=a(1)+a(3);
% %             a(2)<=0;
% %         cvx_end
% %         target=sum(-w*gama);
% % % tar=sum(exp(-abs(a'*n)/sigma));
% %         cvx_clear
% %         cvx_begin 
% %         cvx_solver sedumi
% %         cvx_precision best
% %         cvx_quiet(1)
% %         variable a1(6,1);
% %         variable gama(N_n,1);
% %         minimize sum(-w*gama)
% %         subject to
% % %             abs(a1'*n)'<=gama;
% %        (n'*a1).^2<=gama;
% %             norm([a1(2)+epsilon;a1(1)-a1(3)])<=a1(1)+a1(3);
% %             a1(2)>=0;
% %         cvx_end
% %         target1=sum(-w*gama);
% %         storetarget(1,it)=target;
% %         storetarget(2,it)=target1;
% % % tar1=sum(exp(-abs(a1'*n)/sigma));
% % % storetar(1,it)=tar;
% % % storetar(2,it)=tar1;
% %             if(target<=target1)
% %                 Ma=a;
% %             else
% %                 Ma=a1;
% %             end
% %         if it ==1
% %             aa = Ma;
% %         else
% %             aa0 = aa;
% %             aa = Ma;
% %             if norm(aa-aa0)<1e-5
% %                 break;
% %             end
% %         end
% %         else
% %             if(aa(2)<0)
% %                 cvx_clear
% %                 cvx_begin 
% %                 cvx_solver sedumi
% %                 cvx_precision best
% %                 cvx_quiet(1)
% %                 variable a(6,1);
% %                 variable gama(N_n,1);
% %                 minimize sum(-w*gama)
% %                 subject to
% % %                     abs(a'*n)'<=gama;
% %                 (n'*a).^2<=gama;
% %                     norm([a(2)-epsilon;a(1)-a(3)])<=a(1)+a(3);
% %                     a(2)<=0;
% %                 cvx_end
% %                 target=sum(-w*gama);
% %                 storetarget(1,it)=target;
% % %         tar=sum(exp(-abs(a'*n)/sigma))
% % %         storetar(1,it)=tar;
% %             else
% %                 cvx_clear
% %                 cvx_begin 
% %                 cvx_solver sedumi
% %                 cvx_precision best
% %                 cvx_quiet(1)
% %                 variable a(6,1);
% %                 variable gama(N_n,1);
% %                 minimize sum(-w*gama)
% %                 subject to
% % %                     abs(a'*n)'<=gama;
% %                 (n'*a).^2<=gama;
% %                     norm([a(2)+epsilon;a(1)-a(3)])<=a(1)+a(3);
% %                     a(2)>=0;
% %                 cvx_end
% %                 target1=sum(-w*gama);
% %                 storetarget(2,it)=target1;
% % %         tar=sum(exp(-abs(a'*n)/sigma))
% % %         storetar(2,it)=tar1;
% %             end
% % 
% %             aa0 = aa;
% %             aa = a;
% %             if norm(aa-aa0)<1e-5
% %                 break;
% %             end
% %         end
% %         aaa2(:,it)=aa;
% %         % % ==============================================================更新权重 
% %         data_points = aa.'*n;
% %         sigma_E = std(data_points);
% % 
% %         X = sort(data_points,'ascend');
% %         Q1=X(floor((N_n)/4));
% %         Q3=X(floor(3*(N_n)/4));
% %         r = (Q3-Q1)/1.34;
% %         sigma = 1.06*min(sigma_E,r)*N_n^(-0.2);
% %         ww = -exp(-(data_points).^2/(2*sigma^2))
% %         w = -ww/sum(ww);
% %     end
% %         aaa2(:,it)=aa;
% %     itLPit=it;
% %     data_points = aa.'*n;
% %     w1=log(1+log(1+log(1+abs(data_points))));
% %     IDX=kmeans(w1',2);
% %     me1=abs(sum(w(IDX==1))/sum(IDX==1));
% %     me2=abs(sum(w(IDX==2))/sum(IDX==2));
% % if(me2<me1)
% %     outlier_number=sum(IDX==2);
% % else
% %     outlier_number=sum(IDX==1);
% % end
% %     if((me2>me1&&sum(IDX==1)<=60)||(me2<me1&&sum(IDX==2)<=60))
% % %        disp('成功');
% %         O=1;
% %     else
% % %         disp('失败');
% %         O=0;
% %     end
% %     if(O==1)
% %         break
% %     end
% %     end 
% % end
end