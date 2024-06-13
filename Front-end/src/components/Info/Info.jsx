import React from 'react';
import PropTypes from 'prop-types';

const Info = ({callable_bond_price, k, maturity_date, delta, pi, coupon_rate, face_value, day_count_convention, id, deleteCard}) => {
  const handleDelete = () => {
    deleteCard(id);
  };

  return (
    <div className="col m6 s12">
      <div className="card">
        <div className="card-content">
          <span className="card-title" data-test="callable_bond_price">
          Callable bond price: {callable_bond_price}
          </span>
          <div className="card-data">
            <span data-test="k">Strike Price: {k} </span>
            <span data-test="maturity_date"> Maturity Date: {maturity_date} </span>
            <span data-test="delta"> Delta: {delta} </span>
            <span data-test="pi"> Pi: {pi} </span>

            <span data-test="coupon_rate"> Coupon Rate: {coupon_rate} </span>
            <span data-test="face_value"> Face Value: {face_value} </span>
            <span data-test="day_count_convention"> Day Count Convention: {day_count_convention.label} </span>
          </div>

          <button className="delete-btn" onClick={handleDelete}>
            X
          </button>
        </div>
      </div>
    </div>
  );
};

Info.propTypes = {
  k: PropTypes.string, 
  maturity_date: PropTypes.string, 
  delta: PropTypes.string, 
  pi: PropTypes.string, 
  coupon_rate: PropTypes.string, 
  face_value: PropTypes.string, 
  day_count_convention: PropTypes.string,
  callable_bond_price: PropTypes.string
};

export default Info;
