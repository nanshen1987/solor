function [in_profile,no_epochs,ok]= Read_profile(filename)
%Read_profile - inputs a motion profile in the following .csv format
% Column 1: time (sec)
% Column 2: latitude (deg)γ��
% Column 3: longitude (deg)����
% Column 4: height (m)
% Column 5: north velocity (m/s)
% Column 6: east velocity (m/s)
% Column 7: down velocity (m/s)
% Column 8: roll angle of body w.r.t NED (deg)�����
% Column 9: pitch angle of body w.r.t NED (deg)������
% Column 10: yaw angle of body w.r.t NED (deg)��ƫ��
%
% Software for use with "Principles of GNSS, Inertial, and Multisensor
% Integrated Navigation Systems," Second Edition.
%
% This function created 31/3/2012 by Paul Groves
%
% Inputs:
%   filename     Name of file to write
%
% Outputs:
%   in_profile   Array of data from the file
%   no_epochs    Number of epochs of data in the file
%   ok           Indicates file has the expected number of columns 
%                �����ļ���������Ԥ�ڵ�һ��

% Copyright 2012, Paul Groves
% License: BSD; see license.txt for details

% Begins

% Parameters
deg_to_rad = 0.01745329252;

% Read in the profile in .csv format
pre_in_profile = dlmread(filename);
[pre_no_epochs,pre_no_columns] = size(pre_in_profile);

in_profile=zeros(pre_no_epochs,10);
in_profile(:,1)=pre_in_profile(:,1)-floor(pre_in_profile(1,1));
in_profile(:,2:3)=pre_in_profile(:,8:9);
in_profile(:,4)=pre_in_profile(:,4);
in_profile(:,5)=pre_in_profile(:,11);
in_profile(:,6)=pre_in_profile(:,10);
in_profile(:,7)=-pre_in_profile(:,12);
in_profile(:,8)=pre_in_profile(:,7);
in_profile(:,9)=pre_in_profile(:,6);
in_profile(:,10)=pre_in_profile(:,5);

% Determine size of file
[no_epochs,no_columns] = size(in_profile);

% Check number of columns is correct (otherwise return)���������Ƿ���ȷ
if no_columns~=10
    disp('Input file has the wrong number of columns')
    ok = false;
else
    ok = true;
    % Convert degrees to radians �ѽǶ�ת��Ϊ����
    in_profile(:,2:3) = deg_to_rad * in_profile(:,2:3);
    in_profile(:,8:10) = deg_to_rad * in_profile(:,8:10);
end %if

% Ends