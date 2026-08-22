import streamlit as st
from streamlit_keplergl import keplergl_static
from keplergl import KeplerGl
import pandas as pd
import json
import time
import pydeck as pdk
pdk.settings.mapbox_api_key = "pk.eyJ1IjoiY2FuZGljZTExMjgiLCJhIjoiY21peTdlaDZsMGF2ZzNlb2JrMG05ajE2dyJ9.0QVdf3-F-3oeNo88BNoa7Q"

# Set title
st.title("民 航 疫 情 监 测 系 统")

# Hide streamlit style
hide_streamlit_style = """
            <style>
            #MainMenu {visibility: hidden;}
            footer {visibility: hidden;}
            </style>
            """
st.markdown(hide_streamlit_style, unsafe_allow_html=True)

# Tabs for real-time and historical flights
tab1, tab2 = st.tabs(["实时航班", "历史航班"])

# Historical flights tab
with tab2:
    # Load data from data folder
    routes = pd.read_csv('data/flights.v3.3.3.csv')

    col1, col2, col3 = st.columns(3)

    # Select infectious disease
    with col1:
        choice = st.selectbox("选择传染病：", 
                              ["出发地猴痘确诊病例", 
                               "出发地新冠确诊病例"])
        routes["病例"] = routes[choice]

    # Select destination city
    with col2:
        city = st.multiselect("选择航班目的地：", routes["目的地城市"].unique(),
                              default=["Beijing", "Shenzhen", "Dalian"])
        routes = routes[routes["目的地城市"].map(lambda x: x in city)]

    # Select airlines
    with col3:
        airlines = st.multiselect("选择航班：", routes["航空公司名称"].unique(), 
                                  default=routes["航空公司名称"].unique()[:1])
        routes = routes[routes["航空公司名称"].isin(airlines)]

    cases_max = routes["病例"].max()
    cases_min = routes["病例"].min()

    # Select case range
    values = st.slider('病例数范围', 0.0, float(cases_max), (0.1, float(cases_max)))
    routes = routes[routes["病例"].between(values[0], values[1])]

    # Load config from config folder
    with open("config/config.json", "r") as f:
        config = json.loads(f.read())

    # Create Kepler map
    map_1 = KeplerGl(height=600, width=800, config=config)
    map_1.add_data(data=routes.copy(), name="航班")
    keplergl_static(map_1)

    st.write("航班数: ", len(routes))

# Real-time flights tab
with tab1:
    data = pd.read_parquet("data/data.parquet")
    data["datetime"] = pd.to_datetime(data["unix_time"], unit="s")
    real_time_cases = pd.read_csv("data/realtime-cases.csv")
    
    real_time_cases['病例'] = real_time_cases['出发地猴痘确诊病例']
    real_time_cases = real_time_cases[real_time_cases["病例"] > 0]

    # Scale cases
    max_cases = real_time_cases["病例"].max()
    min_cases = real_time_cases["病例"].min()
    real_time_cases["exits_radius"] = real_time_cases["病例"].map(lambda x: (x - min_cases) / (max_cases - min_cases) * 100 + 1)

    # Create scatter plot layer for real-time cases
    cases_layer = pdk.Layer(
        "ScatterplotLayer",
        real_time_cases,
        pickable=True,
        opacity=0.8,
        stroked=True,
        filled=True,
        radius_scale=1500,
        radius_min_pixels=1,
        radius_max_pixels=300,
        line_width_min_pixels=1,
        get_position=["src_lng", "src_lat"],
        get_radius="exits_radius",
        get_fill_color=[255, 0, 0],
        get_line_color=[0, 0, 0],
    )

    # Real-time data loop
    t = data.datetime.drop_duplicates().sort_values()
    placeholder = st.empty()
    for i in range(0, len(t)):
        states = data[data["datetime"] == t.iloc[i]]
        states = states.fillna(0)
        states["roll"] = 0
        states["yaw"] = -states["true_track"]
        states["pitch"] = 90

        color_map_by_country = {
            "Guatemala": (0, 255, 0, 255),
            "Slovenia": (255, 128, 0, 255),
            "Norway": (255, 0, 255, 255),
            "United States": (255, 0, 0, 255),
        }

        states["color"] = states["origin_country"].map(lambda country: color_map_by_country.get(country, (255, 255, 255, 255)))

        # Scenegraph URL for 3D airplane models
        SCENEGRAPH_URL = "https://raw.githubusercontent.com/visgl/deck.gl-data/master/examples/scenegraph-layer/airplane.glb"

        # Layer for scenegraph
        layer = pdk.Layer(
            type="ScenegraphLayer",
            id="scenegraph-layer",
            data=states,
            pickable=True,
            scenegraph=SCENEGRAPH_URL,
            get_position=["longitude", "latitude", "geo_altitude"],
            sizeMinPixels=0.1,
            sizeMaxPixels=1.5,
            get_orientation=["roll", "yaw", "pitch"],
            get_color="color",
            size_scale=550,
        )

        # Set view for map
        view = pdk.ViewState(latitude=31, longitude=114, zoom=3.5, pitch=0, bearing=0)

        with placeholder.container():
            # Display Beijing time
            import datetime
            beijing_time = datetime.datetime.utcnow() + datetime.timedelta(hours=8)
            st.write("北京时间：", beijing_time)

            # Render map with real-time and historical data
            r = pdk.Deck(
    layers=[cases_layer, layer],
    initial_view_state=view,
    map_provider="mapbox",
    # 明确指定卫星图样式
    map_style="mapbox://styles/mapbox/satellite-v9",
    # 明确把 token 传给前端
    api_keys={"mapbox": pdk.settings.mapbox_api_key},
    tooltip={"text": "航班来自: {origin_country} \n ICAO代码: {icao24}"}
)

            st.pydeck_chart(r)

        time.sleep(0.1)
    time.sleep(0.5)
