% extract peak heating and cooling demands for each building units, usuful
% for sizing the HVAC systems (W, W/m2)
% displays specific thermal energy requirements (kWh/m2,y)
% plots operative temperature of each thermal zone
function [output,input] = analyze_buildings(path_name,project_name)
    % read the input file
    input=import_json_file(path_name+"//"+project_name+".json");
    % read the output file(s)
    if (input.multi_year==0)
        output(1).data=import_output_file(path_name+"//"+project_name+"_output.json");
        output(1).year=[];
    else
        for i=1:numel(input.years)
            output(i).data=import_output_file(path_name+"//"+project_name+"_output.json"+"."+year(i));
            output(i).year=input.years(i);
        end    
    end 
    fprintf(1,"Nomenclature\n");
    fprintf(1,"-Q       thermal energy\n");
    fprintf(1,"-peak    peak thermal power\n");
    fprintf(1,"-year    yearly thermal energy\n");
    fprintf(1,"-s       specific per unit of useful floor area\n");
    fprintf(1,"-H       space heating + domestic hot water\n");
    fprintf(1,"-C       space cooling\n");
    fprintf(1,"-DHW     domestic hot water\n");
    fprintf(1,"----------------------------------------------------------------------------------------------------------------------\n");
    fprintf(1,"Building Unit Q_peak_H Q_peak_C Q_peak_s_H Q_peak_s_C   Q_year_H   Q_year_C Q_year_s_H Q_year_s_C Q_DHW/Q_H floor_area\n");
    fprintf(1,"#        #          kW       kW       W/m2       W/m2      kWh/y      kWh/y   kWh/m2,y   kWh/m2,y         -         m2\n");
    fprintf(1,"----------------------------------------------------------------------------------------------------------------------\n");
    tot_peak_heating=0;
    tot_peak_cooling = 0;
    tot_year_heating = 0;
    tot_year_cooling = 0;
    tot_peak_specific_heating = 0;
    tot_peak_specific_cooling = 0;
    tot_year_specific_heating = 0;
    tot_year_specific_cooling = 0;
    tot_year_dhw = 0;
    tot_floor_useful_area = 0;
    for b=1:numel(output(1).data.buildings)
        for u=1:numel(output(1).data.buildings(b).units)
            peak_heating = max(+(output(1).data.buildings(b).units{u}.hvac.Q_heat_hyd+...
                                 output(1).data.buildings(b).units{u}.hvac.Q_heat_AHU+...
                                +output(1).data.buildings(b).units{u}.dhw.Q_dhw));
            peak_cooling = max(-(output(1).data.buildings(b).units{u}.hvac.Q_cool_hyd+...
                                 output(1).data.buildings(b).units{u}.hvac.Q_cool_AHU));
            year_dhw=sum(+output(1).data.buildings(b).units{u}.dhw.Q_dhw);
            year_heating = sum(+(output(1).data.buildings(b).units{u}.hvac.Q_heat_hyd+...
                                 output(1).data.buildings(b).units{u}.hvac.Q_heat_AHU+...
                                +output(1).data.buildings(b).units{u}.dhw.Q_dhw));
            year_cooling = sum(-(output(1).data.buildings(b).units{u}.hvac.Q_cool_hyd+...
                                 output(1).data.buildings(b).units{u}.hvac.Q_cool_AHU));
            floor_useful_area=output(1).data.buildings(b).units{u}.geometry.floor_area_useful;
            peak_specific_heating =peak_heating/floor_useful_area;
            peak_specific_cooling =peak_cooling/floor_useful_area;
            year_specific_heating =year_heating/floor_useful_area;
            year_specific_cooling =year_cooling/floor_useful_area;
            tot_peak_heating=tot_peak_heating+peak_heating;
            tot_peak_cooling=tot_peak_cooling+peak_cooling;
            tot_peak_specific_heating = tot_peak_specific_heating + peak_specific_heating*floor_useful_area;
            tot_peak_specific_cooling = tot_peak_specific_cooling + peak_specific_cooling*floor_useful_area;
            tot_year_specific_heating = tot_year_specific_heating + year_specific_heating*floor_useful_area;
            tot_year_specific_cooling = tot_year_specific_cooling + year_specific_cooling*floor_useful_area;
            tot_year_heating = tot_year_heating + year_heating;
            tot_year_cooling = tot_year_cooling + year_cooling;
            tot_year_dhw = tot_year_dhw + year_dhw;
            tot_floor_useful_area = tot_floor_useful_area + floor_useful_area;
            fprintf(1,"%-8s %-4i %8.0f %8.0f %10.0f %10.0f %10.0f %10.0f %10.0f %10.0f %9.3f %10.0f\n",...
                             output(1).data.buildings(b).name,u,...
                             peak_heating*1e-3,peak_cooling*1e-3,...
                             peak_specific_heating,peak_specific_cooling,...
                             year_heating*1e-3,year_cooling*1e-3,...
                             year_specific_heating*1e-3,year_specific_cooling*1e-3,year_dhw/year_heating,floor_useful_area);
            figure;
            plot(output(1).data.buildings(b).units{u}.thermal_zone.T_op);
            title(sprintf("T operative (°C) building %s unit % i",output(1).data.buildings(b).name,u));
        end
    end
    tot_peak_specific_heating=tot_peak_specific_heating/tot_floor_useful_area;
    tot_peak_specific_cooling=tot_peak_specific_cooling/tot_floor_useful_area;
    tot_year_specific_heating=tot_year_specific_heating/tot_floor_useful_area;
    tot_year_specific_cooling=tot_year_specific_cooling/tot_floor_useful_area;
    
    fprintf(1,"----------------------------------------------------------------------------------------------------------------------\n");
    fprintf(1,"%-8s %-4i %8.0f %8.0f %10.0f %10.0f %10.0f %10.0f %10.0f %10.0f %9.3f %10.0f\n",...
                     "Tot",0,...
                     tot_peak_heating*1e-3,tot_peak_cooling*1e-3,...
                     tot_peak_specific_heating,tot_peak_specific_cooling,...
                     tot_year_heating*1e-3,tot_year_cooling*1e-3,...
                     tot_year_specific_heating*1e-3,tot_year_specific_cooling*1e-3,tot_year_dhw/tot_year_heating,tot_floor_useful_area);
    fprintf(1,"----------------------------------------------------------------------------------------------------------------------\n");
end