function [] = plot_qext_data_vs_model(qext_train,qext_model_train,qext_predict,qext_model_predict,data_qext,x_axis_start_ticks,x_axis_label,number_of_plots_per_rows,number_of_plots_per_columns,legend_label)

  %% Get number of training and prediction data
  [n_qext,N_train] = size(qext_train);
  N_predict = size(qext_predict,2); 
  mets_meas_in_qext = data_qext.mets_meas_in_qext;
  N = N_train + N_predict;

  %% Get number of models to be plotted
  n_models_to_be_plotted = size(qext_model_train,3);

  %% Linestyles
  Linestyle_all = ["r-";"k--";"g:";"y-."];

  %% Plotting the results


  % Build the x-axis label
  % Plot the qext data and the qext computed from column generation
  set(0,'DefaultFigureWindowStyle','normal');
  number_subplots = number_of_plots_per_columns*number_of_plots_per_rows;
  for ii= 1:n_qext
    j = floor((ii-1)/number_subplots);    
    figure(1+j);
    set(gcf,'unit','normalized','position',[0.1,0.1,0.5,0.6])
    if(mod(ii,number_subplots) == 1)
      tiledlayout(number_of_plots_per_rows,number_of_plots_per_columns, "TileSpacing", "compact");
    end
    nexttile
    h = plot(1:N_train,qext_train(ii,:),'bo','linewidth',1);
    if(N_predict > 0)
      hold on
      plot((N_train+1):(N_train+N_predict),qext_predict(ii,:),'rs','linewidth',1);
    end    
    for lll = 1:n_models_to_be_plotted      
      if(N_predict > 0)
        hold on  
        plot(1:(N_train+N_predict),[qext_model_train(ii,:,lll),qext_model_predict(ii,:,lll)],Linestyle_all(lll),'linewidth',1)
      else
        hold on  
        plot(1:N_train,qext_model_train(ii,:,lll),Linestyle_all(lll),'linewidth',1)
      end
    end
    grid on
    ax = ancestor(h,'axes');
    x_rule = ax.XAxis;
    y_rule = ax.YAxis;
    ax.XTick = x_axis_start_ticks;
    ax.XTickLabel = x_axis_label;
    x_rule.FontSize = 8;
    y_rule.FontSize = 8;
    xlim([1,N])
    ylim padded
    title(mets_meas_in_qext(ii),'interpreter','latex')
    xL = xlabel('Experiment','interpreter','latex');
    xL.FontSize = 9;
    grid on
    hold off
    if(mod(ii,number_subplots) == 1)
      if(N_predict > 0)
        legend(["Training data";"Prediction data";legend_label],'Location','best','FontSize',10)
      else
        legend(["Training data";legend_label],'Location','best','FontSize',10)
      end
    end
  end

end