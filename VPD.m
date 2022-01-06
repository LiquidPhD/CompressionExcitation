function [vec_phase_diff] = VPD(IQData)
IQData = IQData;
% IQData = IQData(1:lowerBound,:,:);
sdl = ones([1 size(IQData,1)]);
VMIQ = permute(IQData,[2 1 3]);
[vector_complex_OCE_data,vec_phase_diff] = ...
    vec_meth_snr(VMIQ,sdl,20);
VMIQ = permute(vector_complex_OCE_data,[2 1 3]);
vec_phase_diff = permute(vec_phase_diff,[2 1 3]);
end