
			addpath('C:\Program Files\MATLAB\R2013b\toolbox\shared\rptgen\@rptgen_ud\@appdata_ud')
getSubclasses('sde','C:\Program Files\MATLAB\R2013b\toolbox\finance\finsupport')
			addpath('C:\Program Files\MATLAB\R2013b\toolbox\shared\rptgen\@rptgen_ud\@appdata_ud')
getSubclasses('sde','C:\Program Files\MATLAB\MATLAB Production Server\R2015a\toolbox\finance')
getSubclasses('plugins.DENSEanalysisPlugin','C:\Users\Neo\Dropbox\denseanalysis-plugin_advanced\+plugins')

supClasses = superclasses('plugins.dense3D_plugin') 
t=plugins.dense3D_plugin;
t=import('plugins.dense3D_plugin.*')
t=meta.package.fromName('tests');
			import plugins.dense3D_plugin.*
			obj = subDENSE3D();
			self.dataObj.computeStrains(obj)
			
			package = meta.package.fromName('DataViewer.DENSE3Dviewer');
			classes = package.ClassList;
			% superclasses('DENSE3DPlugin4CrescentOrgan')
			% meta.class.fromName('plugins.DENSEanalysisPlugin')
			% which('')

			tmp = meta.package.getAllPackages;
			ind = cellfun(@(x)strcmpi(x,'3D DENSE Analysis'),{tmp.Name});
			package = meta.class.fromName('plugins.DENSEanalysisPlugin');
			package = meta.class.fromName('DENSE3DPlugin');
			package = meta.class.fromName('DataViewer');
			package = meta.class.fromName('hgsetget');
			package = meta.class.fromName('DENSE3D');
			package.PropertyList(1, 1).SetAccess = 'public'

			isa(handles.hmanager.Plugins(POI).('datalistener'), 'DENSE3D')
			isa(handles.hmanager.Plugins(POI).('hdense'), 'DENSE3D')
			fields = fieldnames(handles.hmanager.Plugins(POI).Handles);
			fields = fieldnames(handles.hmanager.Plugins(POI));
			ind = cellfun(@(x)isa(handles.hmanager.Plugins(POI).(x), 'DENSE3D'),fields);

findall(self.hfig(1),'type','class')
Hmatch = findobj(self.hfig(1),'-class','DENSE3DPlugin')
Hmatch = findobj(self.hfig(1),'-isa','DENSE3DPlugin');
feval('plugins.dense3D_plugin.DENSE3DPlugin')

				self.dataObj.preview();
				% self.dataObj.parameterize();
				
				ind = zeros(size(faces,1),1); ind(basepoints) = 1;
				hmeshes = findall(allchild(gca),'-property','FaceColor');
				set(hmeshes,'FaceColor', 'w');
				set(hmeshes(2), 'FaceColor', 'flat', 'FaceVertexCData', ind, 'CDataMapping', 'scaled');
				
				hold on; scatter3(vertices(basepoints,1),vertices(basepoints,2),vertices(basepoints,3),10,'b','*');
			uicontrol(...
				'style', 'slider',...
				'Max', 1e3,...
                'Units', 'normalized', ...
				'Position', [0.82 0.95 0.12 0.02],...
                'Value', 0);
				

Warning: The shifted operator has small reciprocal condition estimate: 0.000000
indicating that sigma is near an exact eigenvalue. The algorithm may not converge
unless you try a new value for sigma.
> In eigs>checkInputs/LUfactorAminusSigmaB at 995
  In eigs>checkInputs at 806
  In eigs at 93
  In +DENSE3D_Plugin_4CrescentOrgan\private\biharmonic_embedding at 42
  In +DENSE3D_Plugin_4CrescentOrgan\private\longitudinalParameterization at 23
  In DENSE3DPlugin4CrescentOrgan>DENSE3DPlugin4CrescentOrgan.pickSlice at 1045
  In DENSE3DPlugin4CrescentOrgan>@(s,e)pickSlice(self)
sprank(U)%rank

Warning: Matrix is close to singular or badly scaled. Results may be inaccurate.
RCOND =  5.551115e-17. 
> In eigs>AminusSigmaBsolve at 1210
  In eigs>@(v)AminusSigmaBsolve(Bmtimes(v)) at 149
  In eigs at 261
  In +DENSE3D_Plugin_4CrescentOrgan\private\biharmonic_embedding at 42
  In +DENSE3D_Plugin_4CrescentOrgan\private\longitudinalParameterization at 23
  In DENSE3DPlugin4CrescentOrgan>DENSE3DPlugin4CrescentOrgan.pickSlice at 1045
  In DENSE3DPlugin4CrescentOrgan>@(s,e)pickSlice(self)
4.594095772844915e-08
eps=2.220446049250313e-16
max(u)
min(u)
v=pp*(dgAsB\u);
max(v)
min(v)
v=L\v;
max(max(U))
min(min(U))
v=U\v;
v(v<1e-7) = eps;

3,8->122
	B = biharmonic_embedding(vertices, faces, 3, 6);
    D = sum(bsxfun(@minus, B(basepoints,:).', permute(B, [2 3 1])).^2, 1);
    dist2base = min(D, [], 2);
				% indices = dist2base(:) < 1e-3;
				% indices = dist2apex < 1e-3;
				indices = D < 1e-3;
				indices = abs(D-.01) < .01;
				indices = abs(D-1) < 1e-3;
				indices = D<.02;
				D(apexIndex)

				idxEndo = 1;
				zeros(0,2)
				import plugins.DENSE3D_Plugin_4CrescentOrgan.*
				cdata = RVlongitudinalParameterization(self.dataObj.EndocardialMesh(idxEndo).vertices, self.dataObj.EndocardialMesh(idxEndo).faces, self.dataObj.Apex(idxEndo,:));
				ans(end+1,1)=numel(find(abs(cdata-.99) < .01));
				ans(end+1,2)=numel(find(abs(cdata-.01) < .01));
				
				set(self.hShowMesh.hmeshes(idxP2M), 'FaceColor', 'interp', 'FaceVertexCData', cdata, 'CDataMapping', 'scaled');
				set(self.hShowMesh.hmeshes(idxP2M), 'FaceColor', 'flat', 'FaceVertexCData', cdata, 'CDataMapping', 'scaled');

				
				findall(allchild(2),'-property','style');
				set(ans,'Value',get(ans,'Value')+1);
				
				fields = {'Longitudinal','Circumferential'}
				for k = 1:numel(fields)
				indices = abs(self.dataObj.Parameterization(2).(fields{k})-.999) < 1e-3;
				numel(find(indices))/numel(self.dataObj.Parameterization(2).(fields{k}))*100
				indices = abs(self.dataObj.Parameterization(2).(fields{k})-.999) < 1e-3;
				numel(find(indices))/numel(self.dataObj.Parameterization(2).(fields{k}))*100
				end
				
				lParam = strains.Parameterization.Longitudinal(segments{end});
				lParam = strains.Parameterization.Longitudinal([segments{1:6}]);
			
				indices = abs(self.dataObj.Parameterization(2).Circumferential-.999) < 1e-3;
				numel(find(indices))/numel(self.dataObj.Parameterization(2).Circumferential)*100
				indices = abs(self.dataObj.Parameterization(2).Circumferential-.999) < 1e-3;
				numel(find(indices))/numel(self.dataObj.Parameterization(2).Circumferential)*100


				
				indices = abs(self.dataObj.Parameterization(2).Longitudinal-.999) < 1e-3;
				numel(find(indices))/numel(self.dataObj.Parameterization(2).Longitudinal)*100

				indices = abs(self.dataObj.Parameterization(2).Longitudinal-.999) < 1e-3;
				numel(find(indices))/numel(self.dataObj.Parameterization(2).Longitudinal)*100

					self.hShowMesh.hROIs(6,1) = plot3(self.hShowMesh.ax,[],[],[],'DisplayName', ' Contours','LineWidth',1,'color', 'y');
				