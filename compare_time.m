function [Q1, Q2, time1, time2] = compare_time(rel_err,trials)
% Performance test comparison between the simultaneous Gauss-Kronrod
% quadrature (quadgk_sim.m) and the native ArrayValued numerical 
% integration quadrature of matlab (integral.m).
%
% input: 
% rel_err: relative error for the calculations
% trials: number of trials for the comparison
% Q1: integration results by integral.m
% Q2: integration results by quadgk_sim.m
% time1: average time spent by integral.m
% time2: average time spent by quadgk_sim.m
% Usage example: 
% [Q1, Q2, t1, t2]=compare_time(1e-14,1000);
% This will integrate the functions given in "test_func.m" 1000 
% times with the relative error of 1e-14 along the triangular 
% integration contour [0]-->[1+1i]-->[1-1i]-->[0] on the complex 
% z-plane and compare the average time requirements and the 
% integration results.


time1=zeros(1,trials);
time2=zeros(1,trials);
for ii=1:trials
%     t=cputime;
    tic;
    Q1=integral(@(z)test_func(z),0,1+1i,'ArrayValued',true,'RelTol',rel_err);
    Q1=Q1+integral(@(z)test_func(z),1+1i,1-1i,'ArrayValued',true,'RelTol',rel_err);
    Q1=Q1+integral(@(z)test_func(z),1-1i,0,'ArrayValued',true,'RelTol',rel_err);
    time1(ii)=toc;
end
for ii=1:trials
    tic;
    Q2=quadgk_sim(@(z)test_func(z),[0,1+1i],rel_err,100,7);
    Q2=Q2+quadgk_sim(@(z)test_func(z),[1+1i,1-1i],rel_err,100,7);
    Q2=Q2+quadgk_sim(@(z)test_func(z),[1-1i,0],rel_err,100,7);
    time2(ii)=toc;
end
figure;
plot(time1,'r');
hold on;
plot(time2,'k');

xlabel('Trial no.');
ylabel('Time (s.)');
time1=mean(time1); 
line([1,trials],[time1 time1],'Color','red','LineStyle','--');
time2=mean(time2);
line([1,trials],[time2 time2],'Color','black','LineStyle','--');
legend('Times for ArrayValued Matlab integral function','Times for quadgk_sim','Average time for ArrayValued Matlab integral function','Avergae time for quadgk_sim');
fprintf('quadgk_sim is x%6.3f faster in %1.0i trials.\nAverage time for quadgk_sim is %3.8f sec.\nAverage time for ArrayValued integral is %3.8fsec.\n',time1/time2,trials,time2,time1);
fprintf('Maximum deviation between the results is: %.8e\n',max(abs(abs(Q1-Q2)./Q1)));
end
