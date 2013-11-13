function mutation_frequency_plot2( input , name )
% conjoint distribution of two binomial distribution of proportion test
% input M * 4
%     M total number of mutations
%     For each mutation, we have four numbers [x1 N1 x2 N2]:
% x1: number of mutations in type1
% N1: total number of cases in type1
% x2: number of mutations in type2
% N2: total number of cases in type2

if nargin == 0
    disp('No input! Use CLL data from Hossein as an example.')
    input=[16    58     18    174
        14    58     22    174
        11    58     12    174
        6    58      2    174    ];
    name = {'TP53','NOTCH1','SF3B1','FAT1'};
    M = 4;
elseif nargin == 1
    [ M tmp_e4 ] = size( input );
    if tmp_e4 ~= 4 || any( input( : , 2 ) < input( : , 1 ) ) || any( input( : , 4 ) < input( : , 3 ) )
        error( 'Input error!' )
    else
        for i  = 1 : M
            name{ i } = ['gene' num2str( i )];
        end
    end
elseif nargin > 2
    error( 'Too many input arguments!' )
else
    [ M tmp_e4 ] = size( input );
    if tmp_e4 ~= 4 || any( input( : , 2 ) < input( : , 1 ) ) || any( input( : , 4 ) < input( : , 3 ) ) || M ~= numel( name )
        error( 'Input error!' )
    end
    
end


cutoff1 = 0.68;
cutoff2 = 0.95;
cutoff3 = 0.5;

N_circle = 10;

basiccolor='rgbmcypolasvfn';
map=vivid(M*(N_circle+3),basiccolor(1:M),[0.5 1]); % light color based on rgbm
%colormap(hot);
%colormap(gray)


for m = 1 : M
    x1=input(m,1);
    N1=input(m,2);
    x2=input(m,3);
    N2=input(m,4);
    p1 = x1/N1;
    p2 = x2/N2;
    
    p1tmp = binopdf( 0:N1 , N1 , p1 ); % density of dimension 1
    p2tmp = binopdf( 0:N2 , N2 , p2 ); % density of dimension 2
    PMATRIX = p1tmp' * p2tmp; % joint density (assume two dimension are independent)
    
    all_p = sort(reshape(PMATRIX,(N1+1)*(N2+1),1),'descend');
    cut_density_p1 = all_p( min(find(cumsum(all_p)>=cutoff1)) );
    cut_density_p2 = all_p( min(find(cumsum(all_p)>=cutoff2)) );
    cut_density_p3 = all_p( min(find(cumsum(all_p)>=cutoff3)) );
    
    cut_density_p = [cut_density_p3 cut_density_p1 cut_density_p2];
    
    max_density_p = all_p(1);
    
    x=repmat( (0:N2)/N2 , N1+1 , 1 );
    y=repmat( ([0:N1]')/N1 , 1 , N2+1 );
    
    [C,h(m)] = contour(x,y,PMATRIX,'LineWidth',1);
    
    hold on
    
     set( h(m) , 'ShowText','off' ,...
         'LevelList',cut_density_p(3):(max_density_p-cut_density_p(3))/N_circle:max_density_p)
     
     Cld = get(h(m), 'Children');
     for j=1:length(Cld)
         if strcmp(get(Cld(j), 'Type'), 'patch')
% %            disp(get(Cld(j),'cdata'))
% %            set(Cld(j),'cdata',((m-1)*(N_circle+1)+j)*ones(size(get(Cld(j),'cdata')))); 
             set(Cld(j),'EdgeColor',map((m-1)*(N_circle+3)+j,:));
%              if j > 1
%                  set(Cld(j),'LineStyle','-','LineWidth',2);
%              end
         end
     end
    
end


%alpha(0.7)

line('LineStyle','--','linewidth',4,'color',[0.5 0.5 0.5])

set( gca , 'xlim' , [0,0.5] ,'ylim' , [0 ,0.5] ,'fontsize',16)


for m = 1 : M
    x1=input(m,1);
    N1=input(m,2);
    x2=input(m,3);
    N2=input(m,4);
    plot( x2/N2 , x1/N1 , 'x' , 'color' , map((m-1)*(N_circle+3)+1,:) ,'linewidth',3)
    %   text( x2/N2 , x1/N1 , name{m} ,'fontsize',16 )
    text(x2/N2 , x1/N1,[' <--------------- ',name{m}],'FontName','Arial','FontSize',18,'color',map((m-1)*(N_circle+3)+1,:))
end