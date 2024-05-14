function [NWindow,trayon_starpoint] = EventResonse(trayontime,cellTS,Li,Inti,binwidth)

if ~exist('Li')
        Li=-3; end

if ~exist('Inti')
       Inti=8; end

if ~exist('binwidth')
       binwidth=0.2; end 
                        LeftBound=Li; RightBound=Li+Inti;  CurveN=1; trayon_starpoint=1;
                        clear NWindow
                        for i=1:length(trayontime) % nosepokeon trayontime
                            timestemp=trayontime(i);
                            if timestemp<-LeftBound
                                trayon_starpoint=trayon_starpoint+1;
                                if i==length(trayontime)
                                    NWindow=[];
                                end
                                continue
                            end
                            WindowCellTS=cellTS(cellTS>(timestemp+LeftBound)&cellTS<(timestemp+RightBound));
                            WindowCellTS=WindowCellTS-timestemp-LeftBound;
                            TWindow=0:binwidth:(RightBound-LeftBound);
                            [WindowCurve,sedges]=histcounts(WindowCellTS,TWindow);
                            NWindow(:,CurveN)=hdFlatWindowSmoothing(WindowCurve,4);
                            CurveN=CurveN+1;

                            %
                            %                     aaa=corrcoef(NWindow);
                            %                      a=triu(aaa);
                            %                      a(isnan(a))=0;
                            %                      a(a==1)=0;
                            %                      CorrInWindow(Li+21,Inti-4)=sum(a(:))/nnz(a);
                        end

end

