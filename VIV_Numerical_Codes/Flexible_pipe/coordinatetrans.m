function y = coordinatetrans()
%COORDINATETRANS Summary of this function goes here
%   第i个单元上两个节点i，i+1的坐标变换
t=zeros(12:12);
for cc=1:12
    t(cc,cc)=1;
end
y=t;
end

