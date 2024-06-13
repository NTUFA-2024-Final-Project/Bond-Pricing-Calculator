import React, { useState, useEffect } from 'react';
import { v4 as uuidv4 } from 'uuid';
import 'materialize-css/dist/css/materialize.min.css';
import './App.css';
import BondForm from '../BondForm/BondForm';
import Info from '../Info/Info';
import { getData, storeData } from '../../helpers/localStorage';
import axios from 'axios';

const App = () => {
  const initialState = () => getData('data') || [];
  const [state, setState] = useState(initialState);
  const [data, setData] = useState({});

  useEffect(() => {
    storeData('data', state);
    const k = state.map(obj => obj.k);
    const maturity_date = state.map(obj => obj.maturity_date);
    const delta = state.map(obj => obj.delta);
    const pi = state.map(obj => obj.pi);

    /* -- Callable Bond Information -- */
    const coupon_rate = state.map(obj => obj.coupon_rate);
    const face_value = state.map(obj => obj.face_value);
    const day_count_convention = state.map(obj => obj.day_count_convention);

    const callable_bond_price = state.map(obj => obj.callable_bond_price);

    let newData = { k, maturity_date, delta, pi, coupon_rate, face_value, day_count_convention, callable_bond_price};
    setData(newData);
  }, [state]);

  const handleChange = val => {
    /* -- basic information -- */
    let k = parseFloat(val.k);
    let maturity_date = val.maturity_date;
    let delta = parseFloat(val.delta);
    let pi = parseFloat(val.pi);

    /* -- Callable Bond Information -- */
    let coupon_rate = parseFloat(val.coupon_rate);
    let face_value = parseFloat(val.face_value);
    let day_count_convention = val.day_count_convention;
  
    axios.post('http://localhost:3001/calculate', {
      k, 
      maturity_date, 
      delta, 
      pi, 
      coupon_rate, 
      face_value, 
      day_count_convention
    })
      .then(response => { 
      setState(prevState => {
        	  /* -- basic information -- */
            val.k = k;
            val.maturity_date = maturity_date;
            val.delta = delta;
            val.pi = pi;

            /* -- Callable Bond Information -- */
            val.coupon_rate = coupon_rate;
            val.face_value = face_value;
            val.day_count_convention = day_count_convention;

            val.callable_bond_price = response.data.callable_bond_price;

            val.id = uuidv4();
            let newVal = [...prevState, val];
            return newVal;
        });
    })
    .catch(error => {
        console.error('Error calculating bond price:', error);
    });
  };

  const handleDelete = id => {
    storeData('lastState', state);
    let newState = state.filter(i => {
      return i.id !== id;
    });
    setState(newState);
  };

  const handleUndo = () => {
    setState(getData('lastState'));
  };

  return (
    <div className='container'>
      <div className='row center'>
        <h2 className='white-text'> Bond Price Calculator </h2>
      </div>

      <div className='row center'>
        <h4 className='white-text'> Description </h4>
        <h5 style={{ textAlign: 'left', color: 'white' }}>
          We using Ho and Lee (1986) approach to build term structure trees, and using that to calculate callable bond price. 
          The input variable are:
        </h5>
        <h5 style={{ textAlign: 'left', color: 'white' }}>
               。 K: Strike Price  
          <br/>。 Maturity Date
          <br/>。 Delta
          <br/>。 Pi
          <br/>。 Callable bond information:
          <br/><span style={{ paddingLeft: '50px' }}>Coupon Rate</span>
          <br/><span style={{ paddingLeft: '50px' }}>Face Value</span>
          <br/><span style={{ paddingLeft: '50px' }}>Day Count Convention (e.g. 30/360) </span>
        </h5>
        <h5 style={{ textAlign: 'left', color: 'white' }}>
               。 Zero Coupon Bond Term Structure
          <br/><span style={{ paddingLeft: '50px' }}>Time to Maturity</span>
          <br/><span style={{ paddingLeft: '50px' }}>price </span>
        </h5>
      </div>      

      <div className='row'>
        <div className='col m12 s12'>
          <BondForm change={handleChange} />
          <div>
            <div className='row center'>
              <h4 className='white-text'> Callable Bond Price </h4>
            </div>
            <div className='data-container row'>
              {state.length > 0 ? (
                <>
                  {state.map(info => (
                    <Info
                      key={info.id}
                      id={info.id}

                      k={info.k}
                      maturity_date={info.maturity_date}
                      delta={info.delta}
                      pi={info.pi}

                      /* -- Callable Bond Information -- */
                      coupon_rate={info.coupon_rate}
                      face_value={info.face_value}
                      day_count_convention={info.day_count_convention}

                      callable_bond_price={info.callable_bond_price}
                      
                      deleteCard={handleDelete}
                    />
                  ))}
                </>
              ) : (
                  <div className='center white-text'>No log found</div>
                )}
            </div>
          </div>
          {getData('lastState') !== null ? (
            <div className='center'>
              <button className='calculate-btn' onClick={handleUndo}>
                Undo
              </button>
            </div>
          ) : (
              ''
            )}
        </div>
      </div>
    </div>
  );
};

export default App;
