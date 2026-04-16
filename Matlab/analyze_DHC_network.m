% displays the energy performance of the DHW network 
% - energy center (E000) 
%   - heat supply, COP, PLR for heating
%   - cool supply, EER, PLR for cooling
% - substation (S001)
%    - heat supply and cool supply to the load
%    - heat across the network heat exchangers (HXs)
%    - heat pump COP 
%    - heat pump PLR 
% plots 
% - TES temperature (if present)
% - DHW storage temperature of each substation
% - pipa_A and pipe_B temperature for each substation
function [output,input] = analyze_DHC_network(path_name,project_name)
    % read the input file
    input=import_json_file(path_name+"//"+project_name+".json");
    % read the output file(s)
    if (input.multi_year==0)
        output(1).data=import_output_file(path_name+"//"+project_name+"_output.json");
        output(1).year=[];
    else
        for i=1:numel(input.years)
            output(i).data=import_output_file(path_name+"//"+project_name+"_output.json"+"."+input.years(i));
            output(i).year=input.years(i);
        end    
    end 
    dhc_network=output(1).data.thermal_grids{1,1};
    fprintf(1,"Nomenclature\n");
    fprintf(1,"-mean    seasonal value\n");
    fprintf(1,"-year    yearly thermal energy\n");
    fprintf(1,"-H       space heating + domestic hot water\n");
    fprintf(1,"-C       space cooling\n");
    fprintf(1,"---------------------------------------------------------------------------\n");
    fprintf(1,"Type     Name Q_year_H Q_year_C PLR_mean_H PLR_mean_C COP_mean_H EER_mean_C\n");
    fprintf(1,"#        #         kWh      kWh          -          -          -          -\n");
    fprintf(1,"---------------------------------------------------------------------------\n");
    for i=1:numel(dhc_network.energy_centres)
        ec=dhc_network.energy_centres{i};
        if (ec.type_centre=='C000')
            year_H=sum(ec.WWHPC.Q_heat);
            year_C=sum(ec.WWHPC.Q_cool);
            mean_COP=mean(ec.WWHPC.Q_heat(ec.WWHPC.Q_heat>0)./ec.WWHPC.W_hpc(ec.WWHPC.Q_heat>0));
            PLR_H=mean(ec.WWHPC.PLR(ec.WWHPC.Q_heat>0));
            mean_EER=mean(-ec.WWHPC.Q_cool(ec.WWHPC.Q_cool<0)./ec.WWHPC.W_hpc(ec.WWHPC.Q_cool<0));
            PLR_C=mean(ec.WWHPC.PLR(ec.WWHPC.Q_cool<0));
            fprintf(1,"%-8s %-4s %8.0f %8.0f %10.3f %10.3f %10.3f %10.3f\n",...
                             ec.type_centre,ec.name,...
                             year_H*1e-3,year_C*1e-3,...
                             PLR_H,PLR_C,...
                             mean_COP,mean_EER);
            figure;
            plot(ec.TES.thermal_nodes(2).T);
            title("T TES (°C)");
        end
        if (ec.type_centre=='C002')
            figure;
            plot(ec.bhe_field.thermal_nodes(1).T);
            title("T borehole (°C)");
        end
    end
    fprintf(1,"---------------------------------------------------------------------------\n");
    fprintf(1,"\n");
    fprintf(1,"------------------------------------------------------------------------------------------------\n");
    fprintf(1,"Type     Name Q_year_H Q_aux_HP Q_year_C PLR_mean_H COP_mean_H Q_year_evap Q_aux_nw_H Q_aux_nw_C\n");
    fprintf(1,"#        #         kWh      kWh      kWh          -          -         kWh         kWh       kWh\n");
    fprintf(1,"------------------------------------------------------------------------------------------------\n");
    tot_Q_year_H=0;
    tot_Q_aux_HP=0;
    tot_Q_year_C=0;
    tot_Q_year_evap=0;
    tot_Q_aux_nw_H=0;
    tot_Q_aux_nw_C=0;
    for i=1:numel(dhc_network.substations)
        ss=dhc_network.substations(i);
        if (ss.type_sub=='S001')
            year_H=sum(ss.HP.Q_cond+ss.HP.Q_aux_h);
            tot_Q_year_H=tot_Q_year_H+year_H;
            year_aux_H=sum(ss.HP.Q_aux_h);
            tot_Q_aux_HP=tot_Q_aux_HP+year_aux_H;
            year_C=sum(ss.HX_cooling.thermal_nodes(1).Qin);
            tot_Q_year_C=tot_Q_year_C+year_C;
            mean_COP=mean(ss.HP.Q_cond(ss.HP.Q_cond>0)./ss.HP.W_comp(ss.HP.Q_cond>0));
            PLR_H=mean(ss.HP.PLR(ss.HP.Q_cond>0));
            year_evap=sum(ss.HP.Q_evap);
            tot_Q_year_evap=tot_Q_year_evap+year_evap;
            aux_nw_H=sum(ss.Q_aux_h);
            tot_Q_aux_nw_H=tot_Q_aux_nw_H+aux_nw_H;
            aux_nw_C=sum(ss.Q_aux_c);
            tot_Q_aux_nw_C=tot_Q_aux_nw_C+aux_nw_C;
            fprintf(1,"%-8s %-4s %8.0f %8.0f %8.0f %10.3f %10.3f %11.0f %10.0f %10.0f\n",...
                             ss.type_sub,ss.name,...
                             year_H*1e-3,year_aux_H*1e-3,year_C*1e-3,...
                             PLR_H,mean_COP,...
                             year_evap*1e-3,aux_nw_H*1e-3,aux_nw_C*1e-3);
            figure;
            plot(ss.DHW_tank.T_storage);
            title(ss.name+" T DHW storage (°C)");
            figure;
            plot(ss.HX_heat_pump.thermal_nodes(1).T);
            title(ss.name+" T hot HX (°C)");
            figure;
            plot(ss.HX_cooling.thermal_nodes(1).T);
            title(ss.name+" T cold HX (°C)");
        end
    end
    fprintf(1,"------------------------------------------------------------------------------------------------\n");
    fprintf(1,"%-8s %-4s %8.0f %8.0f %8.0f %10.3f %10.3f %11.0f %10.0f %10.0f\n",...
                "","Tot.",tot_Q_year_H*1e-3,tot_Q_aux_HP*1e-3,tot_Q_year_C*1e-3,0,0,tot_Q_year_evap*1e-3,tot_Q_aux_nw_H*1e-3,tot_Q_aux_nw_C*1e-3);
    fprintf(1,"------------------------------------------------------------------------------------------------\n");

    % trend of borehole temperatures
    if (input.multi_year==1 && strcmp(ec.type_centre,'C002'))
        figure;
        hold on
        for i=1:numel(input.years)
            dhc_network=output(i).data.thermal_grids{1,1};
            ec=dhc_network.energy_centres{1};
            plot(ec.bhe_field.thermal_nodes(1).T);
            title("T borehole (°C)");
        end    
        hold off
        legend_labels = arrayfun(@num2str, input.years, 'UniformOutput', false);
        legend(legend_labels);        
    end     

    % trend of substation 1
    if (input.multi_year==1 && strcmp(ss.type_sub,'S001'))
        for i=1:numel(input.years)
            dhc_network=output(i).data.thermal_grids{1,1};
            ss=dhc_network.substations(1);
            v_year_H(i)=sum(ss.HP.Q_cond+ss.HP.Q_aux_h);
            v_year_aux_H(i)=sum(ss.HP.Q_aux_h);
            v_year_C(i)=sum(ss.HX_cooling.thermal_nodes(1).Qin);
            v_mean_COP(i)=mean(ss.HP.Q_cond(ss.HP.Q_cond>0)./ss.HP.W_comp(ss.HP.Q_cond>0));
            v_PLR_H(i)=mean(ss.HP.PLR(ss.HP.Q_cond>0));
            v_year_evap(i)=sum(ss.HP.Q_evap);
            v_aux_nw_H(i)=sum(ss.Q_aux_h);
            v_aux_nw_C(i)=sum(ss.Q_aux_c);
        end    
        figure;
        bar(input.years, 1e-3*(v_year_H+v_year_aux_H))
        title("Substation 1 - heating (kWh/year)")
        xlabel('Year');        
        figure;
        bar(input.years, v_mean_COP)
        title("Substation 1 - HP COP")
        xlabel('Year');        
        figure;
        bar(input.years, 1e-3*v_year_C)
        title("Substation 1 - cooling (kWh/year)")
        xlabel('Year');        
    end     
    
end