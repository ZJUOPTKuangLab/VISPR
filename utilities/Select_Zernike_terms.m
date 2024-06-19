function aberrations = Select_Zernike_terms(zernikenum, standard) % standard:Wyant, OSA

if strcmp(standard,'Wyant')
    switch(zernikenum)
        case 21
            aberrations = [
            % 0	0	0
            % 1	1	0
            % 1	-1	0
            % 2	0	0
            2	2	0   % VAst
            2	-2	0   % OAst
            3	1	0   % Coma x
            3	-1	0   % Coma y
            4	0	0   % Primary sph
            3	3	0   % Trefoil x
            3	-3	0   % Trefoil y
            4	2	0   % Secondary VAst
            4	-2	0   % Secondary OAst
            5	1	0   % Secondary coma x
            5	-1	0   % Secondary coma y
            6	0	0   % Secondary sph
            4	4	0   % Tetrafoil x
            4	-4	0   % Tetrafoil y
            5	3	0   % Secondary trefoil x
            5	-3	0   % Secondary trefoil y
            6	2	0   % Tertiary VAst
            6	-2	0   % Tertiary OAst
            7	1	0   % Tertiary coma x
            7	-1	0   % Tertiary coma y
            8	0	0]; % Tertiary sph
    
        case 32
         aberrations = [
            % 0	0	0
            % 1	1	0
            % 1	-1	0
            % 2	0	0
            2	2	0
            2	-2	0
            3	1	0
            3	-1	0
            4	0	0
            3	3	0
            3	-3	0
            4	2	0
            4	-2	0
            5	1	0
            5	-1	0
            6	0	0
            4	4	0
            4	-4	0
            5	3	0
            5	-3	0
            6	2	0
            6	-2	0
            7	1	0
            7	-1	0
            8	0	0
            5	5	0
            5	-5	0
            6	4	0
            6	-4	0
            7	3	0
            7	-3	0
            8	2	0
            8	-2	0
            9	1	0
            9	-1	0
            10	0	0];
        case 45
            aberrations = [
            % 0	0	0
            % 1	1	0
            % 1	-1	0
            % 2	0	0
            2	2	0
            2	-2	0
            3	1	0
            3	-1	0
            4	0	0
            3	3	0
            3	-3	0
            4	2	0
            4	-2	0
            5	1	0
            5	-1	0
            6	0	0
            4	4	0
            4	-4	0
            5	3	0
            5	-3	0
            6	2	0
            6	-2	0
            7	1	0
            7	-1	0
            8	0	0
            5	5	0
            5	-5	0
            6	4	0
            6	-4	0
            7	3	0
            7	-3	0
            8	2	0
            8	-2	0
            9	1	0
            9	-1	0
            10	0	0
            6	6	0
            6	-6	0
            7	5	0
            7	-5	0
            8	4	0
            8	-4	0
            9	3	0
            9	-3	0
            10	2	0
            10	-2	0
            11	1	0
            11	-1	0
            12	0	0];
        case 60 
            aberrations = [
            % 0	0	0
            % 1	1	0
            % 1	-1	0
            % 2	0	0
            2	2	0
            2	-2	0
            3	1	0
            3	-1	0
            4	0	0
            3	3	0
            3	-3	0
            4	2	0
            4	-2	0
            5	1	0
            5	-1	0
            6	0	0
            4	4	0
            4	-4	0
            5	3	0
            5	-3	0
            6	2	0
            6	-2	0
            7	1	0
            7	-1	0
            8	0	0
            5	5	0
            5	-5	0
            6	4	0
            6	-4	0
            7	3	0
            7	-3	0
            8	2	0
            8	-2	0
            9	1	0
            9	-1	0
            10	0	0
            6	6	0
            6	-6	0
            7	5	0
            7	-5	0
            8	4	0
            8	-4	0
            9	3	0
            9	-3	0
            10	2	0
            10	-2	0
            11	1	0
            11	-1	0
            12	0	0
            7	7	0
            7	-7	0
            8	6	0
            8	-6	0
            9	5	0
            9	-5	0
            10	4	0
            10	-4	0
            11	3	0
            11	-3	0
            12	2	0
            12	-2	0
            13	1	0
            13	-1	0
            14	0	0];
    end
elseif strcmp(standard,'OSA')
     switch(zernikenum)
          case 17
            aberrations = [
                2,  -2,     0; % vertical astigmatism 
                2,  2,    0;  % oblique astigmatism

                3, -3,    0; % vertical trefoil
                3, -1,    0;  % vertical coma
                3,  1,    0;  % Horizontal coma
                3,  3,    0; % Oblique trefoil
        
                4,  -4,   0;
                4,  -2,   0;
                4,   0,   0;    % primary spherical
                4,  2,    0;
                4,  4,    0;
        
                5,  -5,   0; 
                5,  -3,   0;
                5,  -1,   0; 
                5,  1,    0; 
                5,  3,    0;
                5,  5,    0;];
         case 24
            aberrations = [
                2,  2,     0; % vertical astigmatism 
                2,  -2,    0;  % oblique astigmatism

                3, -3,    0; % vertical trefoil
                3, -1,    0;  % vertical coma
                3,  1,    0;  % Horizontal coma
                3,  3,    0; % Oblique trefoil
        
                4,  -4,   0;
                4,  -2,   0;
                4,   0,   0;    % primary spherical
                4,  2,    0;
                4,  4,    0;
        
                5,  -5,   0; 
                5,  -3,   0;
                5,  -1,   0; 
                5,  1,    0; 
                5,  3,    0;
                5,  5,    0;
        
                6,  -6,   0.0;
                6,  -4,   0;
                6,  -2,   0.0; 
                6,   0,   0;  
                6,  2,    0.0;
                6,  4,    0.0;
                6,  6,    0.0;];
         case 41
            aberrations = [
                    2,  2,     0; % vertical astigmatism 
                    2,  -2,    0;  % oblique astigmatism
            
                    3, -3,    0; % vertical trefoil
                    3, -1,    0;  % vertical coma
                    3,  1,    0;  % Horizontal coma
                    3,  3,    0; % Oblique trefoil
            
                    4,  -4,   0;
                    4,  -2,   0;
                    4,   0,   0;    % primary spherical
                    4,  2,    0;
                    4,  4,    0;
            
                    5,  -5,   0; 
                    5,  -3,   0;
                    5,  -1,   0; 
                    5,  1,    0; 
                    5,  3,    0;
                    5,  5,    0;
            
                    6,  -6,   0.0;
                    6,  -4,   0;
                    6,  -2,   0.0; 
                    6,   0,   0;  
                    6,  2,    0.0;
                    6,  4,    0.0;
                    6,  6,    0.0;
            
                    7,  -7, 0.0;
                    7,  -5, 0.0; 
                    7,  -3, 0.0;
                    7,  -1, 0.0; 
                    7,  1,  0.0; 
                    7,  3,  0.0;
                    7,  5,  0.0;
                    7,  7,  0.0;
            
                    8,  -8, 0.0;
                    8,  -6, 0.0;
                    8,  -4, 0.0;
                    8,  -2, 0.0; 
                    8,   0,   0;
                    8,  2,  0.0;
                    8,  4,  0.0;
                    8,  6,  0.0;
                    8,  8,  0.0;];
         
            case 62
                aberrations = [
                        2,  2,     0; % vertical astigmatism 
                        2,  -2,    0;  % oblique astigmatism
                
                        3, -3,    0; % vertical trefoil
                        3, -1,    0;  % vertical coma
                        3,  1,    0;  % Horizontal coma
                        3,  3,    0; % Oblique trefoil
                
                        4,  -4,   0;
                        4,  -2,   0;
                        4,   0,   0;    % primary spherical
                        4,  2,    0;
                        4,  4,    0;
                
                        5,  -5,   0; 
                        5,  -3,   0;
                        5,  -1,   0; 
                        5,  1,    0; 
                        5,  3,    0;
                        5,  5,    0;
                
                        6,  -6,   0.0;
                        6,  -4,   0;
                        6,  -2,   0.0; 
                        6,   0,   0;  
                        6,  2,    0.0;
                        6,  4,    0.0;
                        6,  6,    0.0;
                
                        7,  -7, 0.0;
                        7,  -5, 0.0; 
                        7,  -3, 0.0;
                        7,  -1, 0.0; 
                        7,  1,  0.0; 
                        7,  3,  0.0;
                        7,  5,  0.0;
                        7,  7,  0.0;
                
                        8,  -8, 0.0;
                        8,  -6, 0.0;
                        8,  -4, 0.0;
                        8,  -2, 0.0; 
                        8,   0,   0;
                        8,  2,  0.0;
                        8,  4,  0.0;
                        8,  6,  0.0;
                        8,  8,  0.0;

                        9,  -9, 0.0;
                        9,  -7, 0.0;
                        9,  -5, 0.0; 
                        9,  -3, 0.0;
                        9,  -1, 0.0; 
                        9,  1,  0.0; 
                        9,  3,  0.0;
                        9,  5,  0.0;
                        9,  7,  0.0;
                        9,  9,  0.0;
                        
                        10, -10,    0.0
                        10, -8, 0.0;
                        10, -6, 0.0;
                        10, -4, 0.0;
                        10, -2, 0.0; 
                        10, 2,  0.0;
                        10, 0,  0.0;
                        10, 4,  0.0;
                        10, 6,  0.0;
                        10, 8,  0.0;
                        10, 10, 0.0;];
     end

end
end