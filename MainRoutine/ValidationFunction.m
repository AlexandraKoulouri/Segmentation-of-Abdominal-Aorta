%define sensitivity and specificity
function [sensitivity,specificity,DICE,Mask]=ValidationFunction(BinaryImage,Img,Mask);

if (size(Mask,1)==1)
    Mask= roipoly(Img);
end
P = length(find(Mask == 1)); %total number of positives
N = length(find(Mask == 0)); %total number of negatives

A = BinaryImage-Mask;
%define False and True positives
FP = find(A==1); %indeces of false positive
FN = find(A==-1);%indeces of false negatives
NumFP = length(FP); %number of f. pos.
NumFN= length(FN); %number of f. neg.
NumTP = P - NumFN; %estimation of t. positives 
NumTN =  N-NumFP;   %number of true negatives
sensitivity = NumTP/P;
specificity = NumTN/N;
A_B=length(find(BinaryImage+Mask==2));
B=length(find(BinaryImage==1));
DICE=2*A_B/(P+B);