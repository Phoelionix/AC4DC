import plotly.graph_objs as go
from plotly.subplots import make_subplots
from IPython.display import clear_output, display
import time
import webbrowser

if __name__ == "__main__":
# Create a new subplot
    fig = make_subplots(rows=1, cols=1)

    # Plot a line chart
    trace = go.Scatter(x=[1, 2, 3], y=[4, 5, 6])
    fig.add_trace(trace)



    # Display the plot
    fig.write_html('plot.html')
    webbrowser.open('plot.html')
    time.sleep(0.5)
    # Clear the previous output
    clear_output(wait=True)


    fig.update_traces(x=[1, 1, 1], y=[4, 5, 6])

    fig.write_html('plot.html')
    webbrowser.open('plot.html')




# import plotly.graph_objs as go
# from plotly.subplots import make_subplots
# from IPython.display import clear_output
# import time

# # Create a new subplot
# fig = make_subplots(rows=1, cols=1)

# # Plot a line chart
# trace = go.Scatter(x=[1, 2, 3], y=[4, 5, 6])
# fig.add_trace(trace)

# fig.show()
# time.sleep(1)
# # Clear the previous output
# clear_output(wait=True)


# # Display the plot


# fig.update_traces(x=[1, 2, 1], y=[4, 3, 6])

# # Clear the previous output


# fig.show()