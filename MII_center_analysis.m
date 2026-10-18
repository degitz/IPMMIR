% Number of MRI slices
noslc = 2;

% Number of uCT pixels
nopix = round(noslc*resdiff);

% Location of starting pixel
uctstart = round(size(uct_i_testout,3)/2) - nopix;

% Location of starting slice
mristart = round(size(mri_i_test,3)/2) - noslc;

% Isolate centers
mricent = Norm01(sum(mri_i_test(:,:,mristart:mristart+noslc),3));
uctcent = Norm01(sum(uct_i_testout(:,:,uctstart:uctstart+nopix),3));

MTcent = MII(uctcent, mricent);