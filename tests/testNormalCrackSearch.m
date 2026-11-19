function testNormalCrackSearch
[x,y] = meshgrid(1:100);
R = -100*exp(-((y-55).^2)/2);
path = [(35:65)' 40*ones(31,1)];
[q,G,mask,info] = refineCrackPathNormalMaximum(path,R);
assert(numel(info.offsets)==61 && info.offsets(31)==0);
assert(all(abs(q(:,2)-55)<1e-12));
assert(isequal(q(:,1),path(:,1)));
assert(isequal([G.Nodes.X G.Nodes.Y],q));
assert(all(mask(55,35:65)));
[q,~,~,info] = refineCrackPathNormalMaximum(path,R,struct('radius',10,'density',2));
assert(numel(info.offsets)==41 && all(abs(q(:,2)-50)<1e-12));
[q,~,~,~] = refineCrackPathNormalMaximum(path,ones(100));
assert(isequal(q,path));
R(:) = NaN;
R(35:45,30:70) = 1;
R(44,30:70) = 4;
[q,~,~,~] = refineCrackPathNormalMaximum(path,R);
assert(all(q(:,2)==44));
% Vertical crack searches horizontally, including its endpoint normals.
R = 100*exp(-((x-55).^2)/2);
path = [40*ones(31,1) (35:65)'];
[q,~,~,~] = refineCrackPathNormalMaximum(path,R);
assert(all(q(:,1)==55) && isequal(q(:,2),path(:,2)));
fprintf('Normal maxima, 61/41 samples, density, endpoints, polarity, ties and ROI tests passed.\n');
end
