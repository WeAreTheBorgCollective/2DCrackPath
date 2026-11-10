function p = phaseFieldParameters(options)
%PHASEFIELDPARAMETERS Parameters for geometry-constrained, mechanical-free PF-DIC.
if nargin<1, options=struct(); end
p=struct('lc',5,'mu',0.1,'grayScale',1,'distanceSmoothing',0.1, ...
    'maxGeometryStep',0.25,'maxBacktracks',12,'outerIterations',3);
fields=fieldnames(options);
for k=1:numel(fields)
    if ~isfield(p,fields{k}), error('phaseFieldParameters:Option','Unknown option: %s',fields{k}); end
    p.(fields{k})=options.(fields{k});
end
for name={'lc','mu','grayScale','distanceSmoothing','maxGeometryStep'}
    validateattributes(p.(name{1}),{'numeric'},{'scalar','real','finite','positive'});
end
for name={'maxBacktracks','outerIterations'}
    validateattributes(p.(name{1}),{'numeric'},{'scalar','real','finite','integer','nonnegative'});
end
end
